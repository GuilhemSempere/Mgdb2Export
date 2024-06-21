package com.traviswheeler.libs;

/**
 * This class is pretty much a clone of NINJA's TreeNode
 * Originally taken from https://github.com/TravisWheelerLab/NINJA, distributed under MIT licence
 * Class was edited in order to remove dependency for LogWriter, which we don't need
 *
 * @author SEMPERE
 */

public class BinaryHeap {
    protected int currentSize;
    protected int nextUnusedPos;
    protected int[] freeSlots;
    protected int freeSlotPos;
    protected int[] keys;
    protected int[] val1s;
    protected int[] val2s;
    protected int[] heapArray;
    protected int maxCapacity = Integer.MAX_VALUE;

    public BinaryHeap(int maxCapacity) {
        this.maxCapacity = maxCapacity;
        init(0, 1000); // Default capacity
    }

    public BinaryHeap() {
        init(0, 1000); // Default capacity
    }
    /**
     * Insert into the priority queue.
     * Duplicates are allowed.
     * @param x the item to insert.
     * @return null, signifying that decreaseKey cannot be used.
     * @throws Exception 
     */
    public int insert( int val1, int key) throws Exception {

        if( currentSize  == keys.length )
            doubleArray( );
        
        //stick input values into a slot in the supporting arrays,         
        int freePos;
    	if (freeSlotPos == -1) {
    		freePos = nextUnusedPos++;
    	} else {
    		freePos = freeSlots[freeSlotPos--];
    	}
    	
    	keys[freePos] = key;
        val1s[freePos] = val1;
        
        //then use that slot as the value stored in the PQ 
        int hole = ++currentSize;
        heapArray[0] = freePos;
        
        // Percolate up
        while( key < keys[ heapArray[hole/2] ] ) {
        	heapArray[hole] = heapArray[ hole/2 ];
        	hole /= 2;
        }
        
        heapArray[hole] = freePos;
        return freePos;
    }
    
    public int insert(int val1, int val2, int key) throws Exception {
        if (currentSize == keys.length) doubleArray();

        int freePos = freeSlotPos == -1 ? nextUnusedPos++ : freeSlots[freeSlotPos--];
        keys[freePos] = key;
        val1s[freePos] = val1;
        val2s[freePos] = val2;
        heapArray[0] = freePos;

        int hole = currentSize;
        while (key < keys[heapArray[hole / 2]]) {
            heapArray[hole] = heapArray[hole / 2];
            hole /= 2;
        }
        heapArray[hole] = freePos;
        return freePos;
    }

    public int[] deleteMin() {
        int minPos = heapArray[1];
        int[] min = {val1s[minPos], val2s[minPos]};
        freeSlots[++freeSlotPos] = heapArray[1];
        heapArray[1] = heapArray[currentSize--];
        percolateDown(1);
        return min;
    }

    public boolean isEmpty() {
        return currentSize == 0;
    }

    public int size() {
        return currentSize;
    }

    public void makeEmpty() {
        currentSize = 0;
        nextUnusedPos = 0;
        freeSlotPos = -1;
    }

    protected void init(int startCount, int size) {
        currentSize = startCount;
        nextUnusedPos = startCount;
        freeSlotPos = -1;
        keys = new int[size];
        val1s = new int[size];
        val2s = new int[size];
        heapArray = new int[size + 1];
        freeSlots = new int[size];
    }

    protected void doubleArray() throws Exception {
        if (keys.length >= maxCapacity) throw new Exception("Heap capacity exceeded");
        int newLength = Math.min(keys.length * 2, maxCapacity);
        int[] newFreeSlots = new int[newLength];
        int[] newKeys = new int[newLength];
        int[] newVal1s = new int[newLength];
        int[] newVal2s = new int[newLength];
        int[] newHeapArray = new int[newLength + 1];
        System.arraycopy(freeSlots, 0, newFreeSlots, 0, freeSlotPos + 1);
        System.arraycopy(keys, 0, newKeys, 0, keys.length);
        System.arraycopy(val1s, 0, newVal1s, 0, val1s.length);
        System.arraycopy(val2s, 0, newVal2s, 0, val2s.length);
        System.arraycopy(heapArray, 0, newHeapArray, 0, heapArray.length);
        keys = newKeys;
        val1s = newVal1s;
        val2s = newVal2s;
        heapArray = newHeapArray;
        freeSlots = newFreeSlots;
    }

    private void percolateDown(int hole) {
        int child, pos = heapArray[hole];
        while (hole * 2 <= currentSize) {
            child = hole * 2;
            if (child != currentSize && keys[heapArray[child + 1]] < keys[heapArray[child]]) child++;
            if (keys[heapArray[child]] < keys[pos]) heapArray[hole] = heapArray[child];
            else break;
            hole = child;
        }
        heapArray[hole] = pos;
    }
}