package fr.cirad.mgdb.exporint.tools.nj;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.Stack;

import com.traviswheeler.ninja.TreeNode;

public class NewickTreeRerooter {

    private static final DecimalFormat decimalFormat;

    static {
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
        decimalFormat = new DecimalFormat("0.########", symbols);
    }

    public TreeNode parseNewick(String newick) {
        Stack<TreeNode> stack = new Stack<>();
        TreeNode root = null;
        TreeNode current = null;
        StringBuilder name = new StringBuilder();
        StringBuilder length = new StringBuilder();
        boolean inLength = false;

        for (char c : newick.toCharArray()) {
            switch (c) {
                case '(':
                    if (current != null) stack.push(current);
                    current = new TreeNode();
                    break;
                case ',':
                    if (name.length() > 0 || length.length() > 0) {
                        TreeNode node = new TreeNode(name.toString(), parseLength(length.toString()));
                        if (current.leftChild == null) {
                            current.leftChild = node;
                        } else {
                            current.rightChild = node;
                        }
                        node.parent = current;
                        name = new StringBuilder();
                        length = new StringBuilder();
                        inLength = false;
                    }
                    break;
                case ')':
                    if (name.length() > 0 || length.length() > 0) {
                        TreeNode node = new TreeNode(name.toString(), parseLength(length.toString()));
                        if (current.leftChild == null) {
                            current.leftChild = node;
                        } else {
                            current.rightChild = node;
                        }
                        node.parent = current;
                        name = new StringBuilder();
                        length = new StringBuilder();
                        inLength = false;
                    }
                    if (!stack.isEmpty()) {
                        TreeNode child = current;
                        current = stack.pop();
                        if (current.leftChild == null) {
                            current.leftChild = child;
                        } else {
                            current.rightChild = child;
                        }
                        child.parent = current;
                    } else {
                        root = current;
                    }
                    break;
                case ':':
                    inLength = true;
                    break;
                case ';':
                    break;
                default:
                    if (inLength) {
                        length.append(c);
                    } else {
                        name.append(c);
                    }
                    break;
            }
        }
        return root;
    }

    private double parseLength(String lengthStr) {
        if (lengthStr.isEmpty()) {
            return Float.MAX_VALUE;
        }
        try {
            return decimalFormat.parse(lengthStr).doubleValue();
        } catch (Exception e) {
            throw new RuntimeException("Failed to parse branch length: " + lengthStr, e);
        }
    }

    public String toNewick(TreeNode node) {
        StringBuffer sb = new StringBuffer();
        node.buildTreeString(sb);
        sb.append(";");
        return sb.toString();
    }

    public TreeNode reroot(TreeNode root, String newRootName) {
        TreeNode newRoot = findNode(root, newRootName);
        if (newRoot == null || newRoot == root) return root; // New root not found or already root

        rerootTree(newRoot);
        newRoot.parent = null;
        return newRoot;
    }

    private TreeNode findNode(TreeNode node, String name) {
        if (node == null) return null;
        if (node.name.equals(name)) return node;
        TreeNode found = findNode(node.leftChild, name);
        if (found == null) {
            found = findNode(node.rightChild, name);
        }
        return found;
    }

    private void rerootTree(TreeNode node) {
        TreeNode current = node;
        TreeNode parent = node.parent;
        node.parent = null;
        TreeNode tempChild;

        while (parent != null) {
            if (parent.leftChild == current) {
                parent.leftChild = parent.rightChild;
                parent.rightChild = current;
            } else {
                parent.rightChild = parent.leftChild;
                parent.leftChild = current;
            }

            tempChild = parent.parent;
            parent.parent = current;
            current = parent;
            parent = tempChild;
        }
    }

    public static void main(String[] args) {
        String newick = "(((((((B001:0.00096951,B003:0.00163278):0.00372365,((((((((B006:0.00066719,B010:0.0003587):0.00008681,B009:0.00030377):0.00019547,B015:0.00087985):0.00008955,B011:0.00101049):0.00234667,B012:0.00201756):0.00552908,B013:0.0018554):0.01884335,B007:0.0211952):0.00568844,B014:0):0.00437491):0.00031,B002:0):0.0003724,CX280:0):0.00003846,B005:0.00009622):0.00000147,B008:0.00034287):0,B004:-0.00000113);";
        String newRootName = "B001";

        NewickTreeRerooter rerooter = new NewickTreeRerooter();
        TreeNode root = rerooter.parseNewick(newick);
        if (root != null) {
            TreeNode newRoot = rerooter.reroot(root, newRootName);
            String newNewick = rerooter.toNewick(newRoot);

            System.out.println("Original Newick: " + newick);
            System.out.println("Rerooted Newick: " + newNewick);
        } else {
            System.out.println("Failed to parse the Newick string.");
        }
    }
}