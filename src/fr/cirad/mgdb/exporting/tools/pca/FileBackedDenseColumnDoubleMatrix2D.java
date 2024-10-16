package fr.cirad.mgdb.exporting.tools.pca;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.List;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.DenseColumnDoubleMatrix2D;

public class FileBackedDenseColumnDoubleMatrix2D extends DenseColumnDoubleMatrix2D {
    private File matrixFile;
    private final long bufferSize = 1024 * 1024 * 10; // 10 MB buffer size
    private final List<MappedByteBuffer> buffers = new ArrayList<>();
    private RandomAccessFile file;
    private FileChannel fileChannel;
    private final long totalSize;

    // Constructor
    public FileBackedDenseColumnDoubleMatrix2D(int rows, int columns, File matrixFile) throws IOException {
        super(0, 0);
        this.matrixFile = matrixFile;
        this.totalSize = (long) rows * columns * Double.BYTES;
        initializeFile(false);
    }

    // Constructor for creating a matrix from existing data
    public FileBackedDenseColumnDoubleMatrix2D(double[][] data, File matrixFile) throws IOException {
        super(0, 0);
        this.rows = data.length;
        this.columns = data[0].length;
        this.matrixFile = matrixFile;
        this.totalSize = (long) rows * columns * Double.BYTES;

        // Initialize the file and map memory buffers, then populate with data
        initializeFile(true);
        writeData(data);
        for (MappedByteBuffer buffer : buffers) {
            buffer.force();
        }
    }

    // Initialize the RandomAccessFile and map the buffers
    private void initializeFile(boolean create) throws IOException {
        if (create) {
            if (!matrixFile.exists()) {
                matrixFile.createNewFile();
            }
        }

        // Open the file in read/write mode
        file = new RandomAccessFile(matrixFile, "rw");
        fileChannel = file.getChannel();

        // Map the buffers to the file
        for (long offset = 0; offset < totalSize; offset += bufferSize) {
            long size = Math.min(bufferSize, totalSize - offset);
            MappedByteBuffer buffer = fileChannel.map(FileChannel.MapMode.READ_WRITE, offset, size);
            buffers.add(buffer);
        }
    }

    // Write the 2D array data into the memory-mapped file
    private void writeData(double[][] data) throws IOException {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                setQuick(i, j, data[i][j]);
            }
        }
    }
    
    public double[][] toArray() {
        double[][] array = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                array[i][j] = getQuick(i, j);
            }
        }
        return array;
    }

    // Override getQuick to read from the file-backed memory-mapped buffers
    @Override
    public double getQuick(int row, int column) {
        int index = row * columns + column;
        long offset = index * Double.BYTES;
        int bufferIndex = (int) (offset / bufferSize);
        long bufferOffset = offset % bufferSize;
        return buffers.get(bufferIndex).getDouble((int) bufferOffset);
    }

    // Override setQuick to write to the file-backed memory-mapped buffers
    @Override
    public void setQuick(int row, int column, double value) {
        int index = row * columns + column;
        long offset = index * Double.BYTES;
        int bufferIndex = (int) (offset / bufferSize);
        long bufferOffset = offset % bufferSize;
        buffers.get(bufferIndex).putDouble((int) bufferOffset, value);
    }

    // Close the RandomAccessFile and FileChannel when done
    public void close() throws IOException {
        if (fileChannel != null) {
            fileChannel.close();
        }
        if (file != null) {
            file.close();
        }
    }
    
    @Override
    public DoubleMatrix2D copy() {
    	return this;
//        try {
//            FileBackedDenseColumnDoubleMatrix2D copy = new FileBackedDenseColumnDoubleMatrix2D(rows, columns, new File(matrixFile.getAbsolutePath() + ".copy"));
//            for (int row = 0; row < rows; row++) {
//                for (int col = 0; col < columns; col++) {
//                    copy.setQuick(row, col, getQuick(row, col));
//                }
//            }
//            return copy;
//        } catch (IOException e) {
//            throw new RuntimeException("Failed to create copy", e);
//        }
    }
    
    @Override
    public double[] elements() {
        double[] elements = new double[rows() * columns()];
        int index = 0;
        for (int col = 0; col < columns(); col++) {
            for (int row = 0; row < rows(); row++) {
                elements[index++] = getQuick(row, col);
            }
        }
        return elements;
    }
}