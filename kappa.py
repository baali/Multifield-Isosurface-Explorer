import pyopencl as cl
import numpy as np
from enthought.tvtk.api import tvtk
import sys
import pylab

if __name__ == '__main__':
    '''
    Initialization and getting all sys-args
    format is kappa.py vtk-filename field-name(1) field-name(2)
    '''
    if (len(sys.argv) != 4):
        print '''Missing Parameters.\nUsage is kappa.py vtk-filename field-name(1) field-name(2)'''
        sys.exit(0)
    bins = 100
    reader = tvtk.UnstructuredGridReader()
    # reader.file_name = 'combustion200x200.vtk'
    reader.file_name = sys.argv[1]
    reader.read_all_scalars = 1
    scalars = sys.argv[2:]
    print scalars

    # this is equivalent to ugrid
    op = reader.get_output()
    op.update()
    # to get all cells
    cells = op.get_cells()
    cell_array = cells.to_array()
    # now every eighth element is number which says how many vertices make a cell so, to get all dimension of cells, cell_array[::9]
    # x, y, z values of all points
    coordinates = op.points.to_array()
    pt_data = op.point_data
    # array for scalar fields
    fields = [pt_data.get_array(scalar).to_array() for scalar in scalars ]
    fs_range = np.array(pt_data.get_array(scalars[0]).range)
    kappaFlag = np.ones(1)
    increment = (fs_range[1] - fs_range[0])/bins
    bins = np.zeros(110, dtype='float32')
    binst = np.zeros(110, dtype='float32')

    indices = np.arange(cells.number_of_cells * 9)
    # We have to get rid of each 9th index
    # id_set = indices[indices[(indices % 9 != 0)]]
    id_cells = cell_array[indices % 9 != 0]
    # j = 112*8+112+1
    # id_set = cells.to_array()[j:j+8]
    # gpu side buffers
    num_threads = 100000
    gpu_coord = np.zeros((num_threads*8, 4))
    gpu_fs = np.zeros(num_threads*8)
    gpu_gs = np.zeros(num_threads*8)
    gpu_hs = np.zeros(num_threads*8)
    num = np.ones(1)*num_threads
    # initializing and creating context for OpenCL
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags

    # GPU side buffers
    buf_coords = cl.Buffer(ctx, mf.READ_ONLY , gpu_coord.nbytes)
    buf_fs = cl.Buffer(ctx, mf.READ_ONLY , gpu_fs.nbytes)
    buf_gs = cl.Buffer(ctx, mf.READ_ONLY , gpu_gs.nbytes)
    buf_hs = cl.Buffer(ctx, mf.READ_ONLY , gpu_hs.nbytes)
    buf_kappaFlag = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=kappaFlag)
    buf_range = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=fs_range)
    buf_increment = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=increment)
    buf_binS = cl.Buffer(ctx, mf.WRITE_ONLY, bins.nbytes * num_threads)
    buf_binsT = cl.Buffer(ctx, mf.WRITE_ONLY, binst.nbytes * num_threads)
    buf_bins = cl.Buffer(ctx, mf.READ_WRITE, bins.nbytes)
    buf_binst = cl.Buffer(ctx, mf.READ_WRITE, binst.nbytes)
    buf_num = cl.Buffer(ctx, mf.READ_ONLY, num.nbytes)

    kernel_string = open('part2.cl').read()
    prg = cl.Program(ctx, kernel_string).build()

    cl.enqueue_write_buffer(queue, buf_bins, bins)
    cl.enqueue_write_buffer(queue, buf_binst, binst)
    cl.enqueue_write_buffer(queue, buf_num, num)

    print len(id_cells)

    for i in range(0, len(id_cells), num_threads * 8):
        if i + (num_threads * 8) > len(id_cells):
            break
        gpu_coord[:,:3] = np.copy(coordinates[id_cells[i:i+(num_threads * 8)]])
        gpu_fs = np.copy(fields[0][id_cells[i:i+(num_threads * 8)]])
        gpu_hs = np.copy(fields[0][id_cells[i:i+(num_threads * 8)]])
        gpu_gs = np.copy(fields[1][id_cells[i:i+(num_threads * 8)]])
        cl.enqueue_write_buffer(queue, buf_coords, gpu_coord)
        cl.enqueue_write_buffer(queue, buf_fs, gpu_fs)
        cl.enqueue_write_buffer(queue, buf_gs, gpu_gs)
        cl.enqueue_write_buffer(queue, buf_hs, gpu_hs).wait()
        prg.part2(queue, (num_threads, ), None, buf_coords, buf_fs, buf_gs, \
                    buf_hs, buf_kappaFlag, buf_range, buf_increment, \
                    buf_binS, buf_binsT)
        prg.summ(queue, (110, ), None, buf_binS, buf_binsT, buf_bins, \
                    buf_binst, buf_num)

    print 'Done for loop'
    gpu_coord = np.zeros(((len(id_cells) - i), 4))
    gpu_coord[:,:3] = np.copy(coordinates[id_cells[i:len(id_cells)]])
    gpu_fs = np.copy(fields[0][id_cells[i:len(id_cells)]])
    gpu_hs = np.copy(fields[0][id_cells[i:len(id_cells)]])
    gpu_gs = np.copy(fields[1][id_cells[i:len(id_cells)]])
    cl.enqueue_write_buffer(queue, buf_coords, gpu_coord)
    cl.enqueue_write_buffer(queue, buf_fs, gpu_fs)
    cl.enqueue_write_buffer(queue, buf_gs, gpu_gs)
    cl.enqueue_write_buffer(queue, buf_hs, gpu_hs).wait()
    num =  np.ones(1) * (len(id_cells) - (i+1))
    cl.enqueue_write_buffer(queue, buf_num, num).wait()

    prg.part2(queue, ((len(id_cells) - i), ), None, buf_coords, buf_fs, buf_gs, \
                  buf_hs, buf_kappaFlag, buf_range, buf_increment, \
                  buf_binS, buf_binsT)
    print 'Done part2'
    prg.summ(queue, (110, ), None, buf_binS, buf_binsT, buf_bins, \
                  buf_binst, buf_num)
    print 'Done sum'
    cl.enqueue_read_buffer(queue, buf_bins, bins).wait()
    iso_values = np.arange(fs_range[0], fs_range[1], increment)
    pylab.plot(iso_values, bins[:100])
    pylab.show()
