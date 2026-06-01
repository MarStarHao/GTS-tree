package queries;

import pbtsrindex.Entry;
import pbtsrindex.Node;
import tools.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.PriorityQueue;

/*
* 实现KNN搜索算法
*
* */
//这是另一类范围查询算法 knn查询算法
public class QueryK_Checkpoint {

    static ArrayList<Integer> checkpoints = new ArrayList<>(); //创建chechpoint集合
    private static int nodeAccess;  //节点访问
    private static int checks;  //checks points

    public QueryK_Checkpoint(){}

    //检查点生成
    private static void generateCheckpoints(int delta) {
        int checkpoint = delta - 1; // 这里要参考论文中的公式
        checkpoints.add(checkpoint);

        while (checkpoint < GlobalConf.tsLength) {
            checkpoint += delta; //检查点+detal
            checkpoints.add(checkpoint);
        }
        checkpoints.remove(checkpoints.size()-1);
    }

    //具体的runQueryk算法实现
    /*
    * 坐标，queryTS ，参数K,参数 epsilon, 参数delta ，以及根节点
    * */
    //这是  二元素组
    public static ArrayList< Tuple2< GeolocatedTS, ArrayList<int[]> > > runQueryK(double[] coords,
                                                                              double[] queryTs,
                                                                              int k, //k参数
                                                                              double epsilonTS, int delta,
                                                                              Node root) throws IOException {

        generateCheckpoints(delta);
        //生成一个二维的数组对象，geoTS 以及arraylist 数组对象
        ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results = new ArrayList<>(); //生成一个arrayList数组对象
        nodeAccess = 0;
        runQueryK(coords, queryTs, k, epsilonTS, delta, root, results); // 在下面的函数中 得到一个knn查询结果
        GlobalConf.totalNodeAccess += nodeAccess;
        GlobalConf.totalChecks += checks;
        return results;
    }

   //KNN查询
    public static ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> runQueryK(double[] qcoords, double[] queryTs,
                                                                              int kappa, double epsilonTS,
                                                                              int delta, Node node,
                                                                              ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) throws IOException {

        PriorityQueue<ElementDistance> queue = new PriorityQueue<>(); //设置一个优先队列 主要放置每个元素的距离
        queue.add(new ElementDistance(node, 0)); //增加每个元素的距离到
        nodeAccess = 0;
        while (queue.isEmpty() == false) {

            ElementDistance ed = queue.poll();  //取出元素
            Node n = ed.n;
            if (n instanceof Entry) { //要取出 geotimeseries

                double[] ts = ((GeolocatedTS) ((Entry) n).getEntry()).getRawTimeSeries(); //获得geoTS序列
                ArrayList<int[]> intervals = new ArrayList<>(); //时间序列间隔

                int start, end = 0; //这里是序列的起始点和终止点
                for (int i = 0; i < checkpoints.size(); i++) { //得到检查点的间隔大小
                    //checks++;
                    int checkpoint = checkpoints.get(i); //获得point
                    if (checkpoint <= end)
                        continue;
                    start = checkpoint; //更新到新的起始点和终止点
                    end = checkpoint;
                    if (Math.abs( ts[checkpoint] - queryTs[checkpoint] ) <= epsilonTS) { //小于epsilonts
                        for (int j = checkpoint - 1; j >= 0; j--) { //遍历检查点
                            //checks++; 这里是check 检查点和查询点之间的距离关系
                            if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) { //如果ts[j]和查询点j的差值 小于epsilonTS
                                start = j; //这个赋给start
                            } else {
                                break;
                            }
                        }
                        for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) { //到下一个checkpoint点，满足ts-query小于eplislon
                            //checks++;
                            if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) { //两者进行差值计算
                                end = j;
                            } else {
                                break;
                            }
                        }

                        if (end - start >= delta)
                            intervals.add(new int[]{start, end});
                    }
                }

                if (intervals.size() > 0)
                    results.add(new Tuple2<>((GeolocatedTS) ((Entry) n).getEntry(), intervals)); //增加这个间隔

                if (results.size() == kappa) //kappa主要干什么？
                    break;

            } else if (n.leaf && !(n instanceof Entry)) { //
                nodeAccess++; //
                long startExecuteTime = System.currentTimeMillis();
                GlobalConf.checks2++;
                ArrayList<GeolocatedTS> contents = n.getData(); //获得数据点
                long totalElapsedExecuteTime = System.currentTimeMillis();
                GlobalConf.tmpCount += totalElapsedExecuteTime - startExecuteTime;

                //取得每个geoTS序列数据
                for (GeolocatedTS geoTS : contents) { //直接获得geotimeseries数据
                    // 下面这行代码是得到了查询点到被查点之间的欧式距离
                    double spatialDistance = Functions.euclideanDistance(geoTS.getCoords(), qcoords); //这是两个位置坐标，分别是geots的坐标和查询序列q的空间坐标
                    ElementDistance eld = new ElementDistance(new Entry<>(geoTS.getCoords(), new double[]{1, 1}, geoTS), spatialDistance);
                    queue.add(eld); //Entry获得了一个空间对象实体 坐标，维度，和geoTS
                }
            } else {
                nodeAccess++; //节点访问的个数
                for (Node c : n.children) {
                    if (!c.leaf) { //非叶子节点的距离计算，但是这个的主要作用是什么啦？
                        double spatialDistance = Functions.distanceNode(c, qcoords); //查询坐标和节点C之间的距离
                        queue.add(new ElementDistance(c, spatialDistance)); //这个主要什么作用
                    } else {
                        int[][][] bitVector = c.getPartBitVector();  //获得叶子结点中的bit向量
                        Tuple2<Integer, Boolean> res = Functions.checkMBTS(queryTs, c, checkpoints, epsilonTS, delta, bitVector); //返回每个点到MBTS边框之间的距离
                        checks += res.x; //这点？
                        if (res.y) {
                            double spatialDistance = Functions.distanceNode(c, qcoords);
                            queue.add(new ElementDistance(c, spatialDistance));
                        }
                    }
                }
            }
        }

        return results;
    }
}