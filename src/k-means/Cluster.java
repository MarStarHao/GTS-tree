package kmeans;

import pbtsrindex.Entry;
import pbtsrindex.Node;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
/*
对应论文中的局部聚类算法，cluster类

在BSTR中使用cluster方法，将不同的序列进行聚类处理
 */
public class Cluster implements Serializable {
    public List<Entry> points;
    public Entry centroid;
    public int id;

    public Cluster(){}

    //Creates a new Cluster
    public Cluster(int id) {
        this.id = id;
        this.points = new ArrayList();
        this.centroid = null;
    }

    public List getPoints() {
        return points;
    }

    public void addPoint(Entry point) {
        points.add(point);
    }

    public Entry getCentroid() {
        return centroid;
    }

    public void setCentroid(Entry centroid) {
        this.centroid = centroid;
    }

    public void clear() {
        points.clear();
    }

    public void plotCluster() {
        System.out.println("[Cluster: " + id + "]");
        System.out.println("[Centroid: " + centroid + "]");
        System.out.println("[Points: \n");
        for (Node p : points) {
            System.out.println(p);
        }
        System.out.println("]");
    }

}
