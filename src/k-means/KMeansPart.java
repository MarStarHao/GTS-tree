package kmeans;

import pbtsrindex.Entry;
import pbtsrindex.Node;
import tools.Functions;
import tools.GeolocatedTS;
import tools.GlobalConf;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
/*
Kmeans算法 论文中的分组代码
不太清楚这部分的主要作用
 */
public class KMeansPart implements Serializable {
    private int numPoints; // 数据点的个数
    private List<Node> points; // 点数
    private List<Cluster> clusters;
    private Node parentNode; //父节点
    private int numClusters;
    private int intervalNo; //间隔的数量

    public KMeansPart(){}

    public KMeansPart(int numClusters, Node node, int i) {
        this.parentNode = node;
        this.numClusters = numClusters;
        this.intervalNo = i; //间隔数量
        this.clusters = new ArrayList();
        if (parentNode.isLeaf()) {
            this.points = node.children;  //判断成孩子节点
        } else {
            this.points = new LinkedList<>(); //
            for (Node child : node.children) {
                for (int j = 0; j < child.partKmeans[i].getNumClusters(); j++) { // 获得聚类的个数
                    //获得了geotimeseries 时间序列 即为在cluster中的中心时间序列条数
                    GeolocatedTS tmp = (GeolocatedTS) child.partKmeans[i].getClusters().get(j).
                            getCentroid().getEntry();
                    GeolocatedTS htsd = new GeolocatedTS(tmp.getRawTimeSeries());
                    Entry point = new Entry(htsd); //计算点实体
                    child.calculateMinMaxPartTS(); //计算最小最大值
                    //point.setPartAbsMinTS(Functions.partitionTS(tmp.getRawTimeSeries()));
                    //point.setPartAbsMaxTS(Functions.partitionTS(tmp.getRawTimeSeries()));
                    point.setPartAbsMinTS(child.getPartAbsMinTS());
                    point.setPartAbsMaxTS(child.getPartAbsMaxTS());
                    points.add(point); //在同一个时刻找到最大值和最小值点
                }
            }
        }
        this.numPoints = points.size();
    }

    //Initializes the process
    public void init() { //初始化簇和簇的中心
        //Create Clusters and set random Centroids
        for (int i = 0; i < numClusters; i++) {
            int random = (int) (Math.random() * this.numPoints); //随机点数
            Cluster cluster = new Cluster(i);
            Node point = points.get(random);
            GeolocatedTS tmp = (GeolocatedTS) ((Entry) point).getEntry();
            GeolocatedTS htsd = new GeolocatedTS(tmp.getRawTimeSeries());
            Entry centroid = new Entry(htsd);
            cluster.setCentroid(centroid); //设定了簇中心
            clusters.add(cluster);
        }
    }

    //The process to calculate the K Means, with iterating method.
    public void calculate() throws CloneNotSupportedException {  //计算k-means
        boolean finish = false;
        int maxIter = 150; //迭代次数

        // Add in new data, one at a time, recalculating centroids with each new one.
        // 每次重新计算一次簇中心
        int count = 0;
        while (!finish) {
            //Clear cluster state
            clearClusters(); //清除所有的簇 相当于再次更新簇中心
            List<Entry> lastCentroids = getCentroids(); //得到所有簇中心

            //Assign points to the closer cluster
            assignCluster(); //重新分配簇中心 这里为什么要进行再次分配节点中心

            //Calculate new centroids.
            calculateCentroids();

            List<Entry> currentCentroids = getCentroids();

            //Calculates total distance between new and old Centroids
            // 重新计算过去簇中心点和新簇的中心点之间距离
            double distance = 0;
            for (int i = 0; i < lastCentroids.size(); i++) {
                //节点之间的距离
                distance += Functions.euclideanDistance(lastCentroids.get(i), currentCentroids.get(i), intervalNo);
            }

            if (distance == 0 || count==maxIter) {
                finish = true;
            }
            count++;
        }
    }

    private void assignCluster() {
       double max = Double.MAX_VALUE;
       // double max = -1;
        double min;
        int cluster = 0;
        double distance;

        for (Node point : points) {
            min = max; // 为啥最大赋值
            //max = -1;
            for (int i = 0; i < numClusters; i++) { //
                Cluster c = clusters.get(i);
                //计算了点到簇中心之间的距离， 但是这里看不出为什么要使用？ 作用是什么？
                // 但是这里实在看不出他们的主要作用是啥？
                distance = Functions.euclideanDistance(point, c.getCentroid(), intervalNo);
                if (distance < min) {
                    min = distance;
                    cluster = i;
                }
                /*
                double[] values1 = ((GeolocatedTS) ((Entry<?>) point).getEntry()).getRawTimeSeries();
                double[] values2 = ((GeolocatedTS) ((Entry<?>) c.getCentroid()).getEntry()).getRawTimeSeries();
                double localScore = 0.0;
                boolean broken = true;
                int start = 0, end = 0;
                for (int j = 0; j < GlobalConf.timePartitions[intervalNo].size(); j++) {

                    if (broken) {
                        start = j;
                    }

                    if (Math.abs(values1[j] - values2[j]) <= GlobalConf.epsilonTSbuild) {
                        localScore++;
                        broken = false;
                        end = j;
                        if (j == GlobalConf.tsLength - 1) {
                            if (end - start >= 2) {
                            //    localScore += 0.5*(end - start);
                            }
                        }
                    } else {
                        broken = true;
                        if (end - start >= 2) {
                            //localScore += 0.5*(end - start);
                        }
                    }
                }

                distance = localScore;
                if (distance > max) {
                    max = distance;
                    cluster = i;
                }*/
            }
            clusters.get(cluster).addPoint((Entry) point); //将这个点增加到clusters中
        }
    }
//计算簇中心位置
//
    private void calculateCentroids() { //计算簇中心
        for (Cluster cluster : clusters) {
            double[] sumDims = new double[GlobalConf.tsLength]; // 变量sumDims其实表示的是时间序列的长度
            List<Entry> list = cluster.getPoints(); // 获得clusters的点的总数
            int nPoints = list.size(); //总数
            if (nPoints == 0)
                continue;

            for (Entry point : list) {
                for (int i = 0; i < GlobalConf.tsLength; i++) {
                    sumDims[i] += ((GeolocatedTS) point.getEntry()).getRawTimeSeries()[i]; //每个点的维度，相当于是序列的长度
                }
            }
            Entry centroid = cluster.getCentroid(); //得到簇中心
            double[] paa = new double[GlobalConf.tsLength];
            for (int i = 0; i < GlobalConf.tsLength; i++) {
                paa[i] = sumDims[i] / nPoints; //相当于每个维度的加和以后取平均，这表示在每个时序点上平均数
            }
            GeolocatedTS htsd = new GeolocatedTS(paa);
            centroid.setEntry(htsd);
        }
    }

    private void clearClusters() {
        for (Cluster cluster : clusters) {
            cluster.clear();
        }
    }
//   获得簇中心
    private List<Entry> getCentroids() throws CloneNotSupportedException {
        List<Entry> centroids = new ArrayList(numClusters);
        for (Cluster cluster : clusters) {
            Entry point = (Entry) cluster.getCentroid().clone();
            centroids.add(point);
        }
        return centroids;
    }

    public int getNumClusters() {
        return numClusters;
    }

    public List<Cluster> getClusters() {
        return clusters;
    }

    public Node getParentNode() {
        return parentNode;
    }
}