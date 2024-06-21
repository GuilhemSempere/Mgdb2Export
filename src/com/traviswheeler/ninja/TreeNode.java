package com.traviswheeler.ninja;

import java.text.DecimalFormatSymbols;
import java.util.Locale;

import java.text.DecimalFormat;

/**
 * This class is pretty much a clone of NINJA's TreeNode
 * A DecimalFormat instance was added in order to force the use of dot as a decimal separator, whatever the JVM locale
 *
 * @author SEMPERE
 */

public class TreeNode {
	public TreeNode parent = null;
	public TreeNode leftChild = null;
	public TreeNode rightChild = null;
	public String name;
	public double length = Float.MAX_VALUE;
	
	private static DecimalFormat df;
	
	static {
		DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.US);
		otherSymbols.setDecimalSeparator('.');
		df = new DecimalFormat("##0.########", otherSymbols);
	}
	
	public TreeNode() {
		name = "";
	}

	public TreeNode(String name) {
		this.name = name;
	}

	public TreeNode(String nodeName, double branchLength) {
		this.name = nodeName;
		this.length = branchLength;
	}

	final public void buildTreeString (StringBuffer sb) {
		String len;
		if (length == Float.MAX_VALUE)
			len = "";
		else 
			len = ":" + df.format(length);
	
		if (null == leftChild) {
			sb.append(name  + len);
		} else {
			sb.append("(");
			leftChild.buildTreeString(sb);
			sb.append(",");
			rightChild.buildTreeString(sb);
			sb.append(")" + len);
		} 
	}	
}