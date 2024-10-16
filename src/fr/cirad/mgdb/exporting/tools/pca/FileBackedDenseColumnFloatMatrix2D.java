package fr.cirad.mgdb.exporting.tools.pca;

import org.ejml.simple.SimpleMatrix;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

public class FileBackedDenseColumnFloatMatrix2D {

    private final int rows;
    private final int columns;
    private MappedByteBuffer buffer;

    public FileBackedDenseColumnFloatMatrix2D(float[][] data, File matrixFile) throws IOException {
        this.rows = data.length;
        this.columns = data[0].length;

        try (RandomAccessFile file = new RandomAccessFile(matrixFile, "rw")) {
            FileChannel channel = file.getChannel();
            long size = rows * columns * Float.BYTES;
            buffer = channel.map(FileChannel.MapMode.READ_WRITE, 0, size);

            // Write the data to the memory-mapped file
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    buffer.putFloat(i * columns * Float.BYTES + j * Float.BYTES, data[i][j]);
                }
            }
        }
    }

    public SimpleMatrix getSimpleMatrix() {
        // Load the data back from the memory-mapped file into a SimpleMatrix (EJML's float-based matrix)
        float[][] data = new float[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                data[i][j] = buffer.getFloat(i * columns * Float.BYTES + j * Float.BYTES);
            }
        }

        return new SimpleMatrix(data);
    }
    
    // Helper method to calculate the correct buffer index
    private int getIndex(int i, int j) {
        return (i * columns + j) * Float.BYTES;
    }

    // Method to get the value at a specific (i, j) position from the memory-mapped file
    public float get(int i, int j) {
        if (i < 0 || i >= rows || j < 0 || j >= columns) {
            throw new IndexOutOfBoundsException("Invalid matrix indices");
        }
        return buffer.getFloat(getIndex(i, j));
    }

	public int numRows() {
		return rows;
	}

	public int numCols() {
		return columns;
	}
}
