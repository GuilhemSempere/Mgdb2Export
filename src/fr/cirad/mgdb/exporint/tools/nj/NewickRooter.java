package fr.cirad.mgdb.exporint.tools.nj;

import java.util.*;
import java.util.regex.*;

public class NewickRooter {
    public static void main(String[] args) {
        String newickTree = "(((((CR1434:0.12357,((CR1899:0.20245,((CR1903:0.01753,(CR2526:0.32592,CR2566:0.38624):0.14991):0.11578,CR2568:0.02011):0.19799):0.084,CR2524:0.16591):0.04444):0.03253,CR2560:0.19256):0.02335,((CR1900:0.21155,CR1904:0.21051):0.02091,CR2577:0.11748):0.07297):0.01279,CR1901:0.18977):0.01257,(CR2528:0.20566,(CR2529:0.23734,CR2569:0.45425):0.03724):0.01257);";
        String wantedRoot = "CR2560";
        
        try {
            String rootedNewickTree = rootTree(newickTree, wantedRoot);
            System.out.println(rootedNewickTree);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static String rootTree(String newickTree, String wantedRoot) {
        // Parse the tree and find the leaf
        Node root = parseNewickTree(newickTree.trim());
        Node targetLeaf = findLeaf(root, wantedRoot);
        
        if (targetLeaf == null) {
            throw new IllegalArgumentException("Specified leaf not found in the tree");
        }
        
        // Root the tree at the target leaf
        List<Node> pathToRoot = getPathToRoot(targetLeaf);
        Node newRoot = reRootTree(pathToRoot);
        
        return newRoot.toNewickString();
    }

    private static Node parseNewickTree(String newickTree) {
        Stack<Node> stack = new Stack<>();
        Node root = null;
        Matcher matcher = Pattern.compile("[(),;]|[^(),;]+").matcher(newickTree);
        
        while (matcher.find()) {
            String token = matcher.group();
            
            switch (token) {
                case "(":
                    stack.push(new Node());
                    break;
                case ")":
                    Node node = stack.pop();
                    if (!stack.isEmpty()) {
                        stack.peek().addChild(node);
                    } else {
                        root = node;
                    }
                    break;
                case ",":
                    break;
                case ";":
                    break;
                default:
                    String[] parts = token.split(":");
                    Node leaf = new Node(parts[0], parts.length > 1 ? Double.parseDouble(parts[1]) : 0.0);
                    if (!stack.isEmpty()) {
                        stack.peek().addChild(leaf);
                    } else {
                        root = leaf;
                    }
                    break;
            }
        }
        
        return root;
    }

    private static Node findLeaf(Node node, String name) {
        if (node.name != null && node.name.equals(name)) {
            return node;
        }
        
        for (Node child : node.children) {
            Node result = findLeaf(child, name);
            if (result != null) {
                return result;
            }
        }
        
        return null;
    }

    private static List<Node> getPathToRoot(Node node) {
        List<Node> path = new ArrayList<>();
        
        while (node != null) {
            path.add(node);
            node = node.parent;
        }
        
        return path;
    }

    private static Node reRootTree(List<Node> pathToRoot) {
        for (int i = 0; i < pathToRoot.size() - 1; i++) {
            Node current = pathToRoot.get(i);
            Node parent = pathToRoot.get(i + 1);
            
            double branchLength = current.branchLength;
            parent.children.remove(current);
            current.children.add(parent);
            current.branchLength = parent.branchLength;
            parent.branchLength = branchLength;
            parent.parent = current;
        }
        
        return pathToRoot.get(0);
    }
}

class Node {
    String name;
    double branchLength;
    List<Node> children;
    Node parent;

    Node() {
        this.children = new ArrayList<>();
    }

    Node(String name, double branchLength) {
        this();
        this.name = name;
        this.branchLength = branchLength;
    }

    void addChild(Node child) {
        child.parent = this;
        children.add(child);
    }

    String toNewickString() {
        if (children.isEmpty()) {
            return name + ":" + branchLength;
        } else {
            StringBuilder sb = new StringBuilder();
            sb.append("(");
            for (int i = 0; i < children.size(); i++) {
                if (i > 0) {
                    sb.append(",");
                }
                sb.append(children.get(i).toNewickString());
            }
            sb.append(")").append(name != null ? name : "").append(":").append(branchLength);
            return sb.toString();
        }
    }
}
