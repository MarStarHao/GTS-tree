package queries;

import pbtsrindex.Node;
import tools.Functions;
import tools.GeolocatedTS;
import tools.GlobalConf;
import tools.Tuple2;

import java.io.IOException;
import java.util.ArrayList;


/*
* range查询，类似于knn查询算法
* */
public class QueryR_Checkpoint {

    static ArrayList<Integer> checkpoints = new ArrayList<>();
    private static int nodeAccess;
    private static int checks;
    private static int spatialPrunes;

    public QueryR_Checkpoint(){}

    private static void generateCheckpoints(int delta) {
        int checkpoint = delta - 1;
        checkpoints.add(checkpoint);

        while (checkpoint < GlobalConf.tsLength) {
            checkpoint += delta;
            checkpoints.add(checkpoint);
        }
        checkpoints.remove(checkpoints.size()-1);
    }
    //------------------------------------------------------
    //𝑄𝑟𝑟(𝑇𝑞, 𝜌, 𝜖, 𝛿)
    public static ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> runQueryR(double[] coords,
                                                                                        double[] queryTs,
                                                                                        double epsilonSP,
                                                                                        double epsilonTS, int delta,
                                                                                        Node root) throws IOException {

        generateCheckpoints(delta);  //生成检查点
        ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results = new ArrayList<>();
        nodeAccess = 0;
        runQueryR(coords, queryTs, epsilonSP, epsilonTS, delta, root, results);
        GlobalConf.totalNodeAccess += nodeAccess;
        GlobalConf.totalChecks += checks;
        GlobalConf.totalSpatialPrunes += spatialPrunes;
        return results;
    }
    //---------------------------------------------------------------------------
    private static void runQueryR(double[] qcoords, double[] queryTs, double epsilonSP, double epsilonTS,
                                  int delta, Node n, ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) throws IOException {

        if (n.leaf) {
            nodeAccess++;            //访问是叶子结点
            long startExecuteTime = System.currentTimeMillis();
            GlobalConf.checks2++;    // 这是什么作用
            ArrayList<GeolocatedTS> contents = n.getData(); //获得节点上时序数据
            long totalElapsedExecuteTime = System.currentTimeMillis();
            GlobalConf.tmpCount += totalElapsedExecuteTime - startExecuteTime;
            for (GeolocatedTS geoTS : contents) {

                double spatialDistance = Functions.euclideanDistance(geoTS.getCoords(), qcoords);
                if (spatialDistance < epsilonSP) {

                    double[] ts = geoTS.getRawTimeSeries();

                    ArrayList<int[]> intervals = new ArrayList<>();

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;
                        if (Math.abs(ts[checkpoint] - queryTs[checkpoint]) <= epsilonTS) {
                            for (int j = checkpoint - 1; j >= 0; j--) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
                                    start = j;
                                } else {
                                    break;
                                }
                            }
                            for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
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
                        results.add(new Tuple2<>(geoTS, intervals));
                }
                else {
                    spatialPrunes++;
                }
            }

        } else {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance > epsilonSP) {
                    spatialPrunes++;
                    continue;
                }

                if (!c.leaf)
                    runQueryR(qcoords, queryTs, epsilonSP, epsilonTS, delta, c, results);
                else {
                    int[][][] bitVector = c.getPartBitVector();
                    Tuple2<Integer, Boolean> res = Functions.checkMBTS(queryTs, c, checkpoints, epsilonTS, delta, bitVector);
                    checks += res.x;
                    if (res.y) {
                        runQueryR(qcoords, queryTs, epsilonSP, epsilonTS, delta, c, results);
                    }
                }
            }
        }
    }

    /*
    private static void runQueryR(double[] qcoords, double[] queryTs, double epsilonSP, double epsilonTS,
                                  int delta, Node n, ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) {
        nodeAccess++;
        if (n.leaf) {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance < epsilonSP) {
                    double[] ts = ((GeolocatedTS) ((Entry) c).getEntry()).getRawTimeSeries();

                    ArrayList<int[]> intervals = new ArrayList<>();

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;
                        if (Math.abs(ts[checkpoint] - queryTs[checkpoint]) <= epsilonTS) {
                            for (int j = checkpoint - 1; j >= 0; j--) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
                                    start = j;
                                } else {
                                    break;
                                }
                            }
                            for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
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
                        results.add(new Tuple2<>((GeolocatedTS) ((Entry) c).getEntry(), intervals));
                }
            }

        } else {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance > epsilonSP)
                    continue;

                boolean found = false;
                for (int j = 0; j < c.getMaxTS().length; j++) {

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;

                        //if ((queryTs[checkpoint] > c.getMaxTS()[j][checkpoint] && queryTs[checkpoint] - c.getMaxTS()[j][checkpoint] <= epsilonTS) ||
                        //        (queryTs[checkpoint] < c.getMinTS()[j][checkpoint] && c.getMinTS()[j][checkpoint] - queryTs[checkpoint] <= epsilonTS) ||
                        //        (queryTs[checkpoint] <= c.getMaxTS()[j][checkpoint] && queryTs[checkpoint] >= c.getMinTS()[j][checkpoint])) {
                        if (verifyPartitions(queryTs[checkpoint], checkpoint, epsilonTS, c)) {

                            for (int k = checkpoint - 1; k >= 0; k--) {
                                checks++;
                                if (queryTs[k] <= c.getMaxTS()[j][k] && queryTs[k] >= c.getMinTS()[j][k]) {
                                    if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                        start = k;
                                    else
                                        break;
                                } else if (queryTs[k] > c.getMaxTS()[j][k]) {
                                    if (queryTs[k] - c.getMaxTS()[j][k] <= epsilonTS) {
                                       // if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                            start = k;
                                      //  else
                                       //     break;
                                    } else {
                                        break;
                                    }
                                } else if (queryTs[k] < c.getMinTS()[j][k]) {
                                    if (c.getMinTS()[j][k] - queryTs[k] <= epsilonTS) {
                                      //  if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                            start = k;
                                      //  else
                                      //      break;
                                    } else {
                                        break;
                                    }
                                }
                            }

                            for (int k = checkpoint + 1; k < GlobalConf.tsLength; k++) {
                                checks++;
                                if (queryTs[k] <= c.getMaxTS()[j][k] && queryTs[k] >= c.getMinTS()[j][k]) {
                                    if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                        end = k;
                                    else
                                        break;
                                } else if (queryTs[k] > c.getMaxTS()[j][k]) {
                                    if (queryTs[k] - c.getMaxTS()[j][k] <= epsilonTS) {
                                       // if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                            end = k;
                                      //  else
                                       //     break;
                                    } else {
                                        break;
                                    }
                                } else if (queryTs[k] < c.getMinTS()[j][k]) {
                                    if (c.getMinTS()[j][k] - queryTs[k] <= epsilonTS) {
                                      //  if (verifyPartitions(queryTs[k], k, epsilonTS, c))
                                            end = k;
                                      //  else
                                      //      break;
                                    } else {
                                        break;
                                    }
                                }
                            }

                            if (end - start >= delta) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (found)
                        break;
                }

                if (found)
                    runQueryR(qcoords, queryTs, epsilonSP, epsilonTS, delta, c, results);

            }
        }
    }

    private static boolean verifyPartitions(double value, int k, double epsilonTS, Node c) {

        int partition = Functions.getPartition(k);
        int bundleIndex = (int) (k%(Math.ceil((double)GlobalConf.tsLength/GlobalConf.numTimePartitions)));
        for (int j = 0; j < GlobalConf.numClusters; j++) {
            checks++;
            if ((value > c.getPartMaxTS()[partition][j][bundleIndex] && value - c.getPartMaxTS()[partition][j][bundleIndex] <= epsilonTS) ||
                    (value < c.getPartMinTS()[partition][j][bundleIndex] && c.getPartMinTS()[partition][j][bundleIndex] - value <= epsilonTS) ||
                    (value <= c.getPartMaxTS()[partition][j][bundleIndex] && value >= c.getPartMinTS()[partition][j][bundleIndex])) {
                return true;
            }
        }

        return false;
    }

    private static void runQueryR(double[] qcoords, double[] queryTs, double epsilonSP, double epsilonTS,
                                            int delta, Node n, ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) {
        nodeAccess++;
        if (n.leaf) {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance < epsilonSP) {
                    double[] ts = ((GeolocatedTS) ((Entry) c).getEntry()).getRawTimeSeries();

                    ArrayList<int[]> intervals = new ArrayList<>();

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;
                        if (Math.abs(ts[checkpoint] - queryTs[checkpoint]) <= epsilonTS) {
                            for (int j = checkpoint - 1; j >= 0; j--) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
                                    start = j;
                                } else {
                                    break;
                                }
                            }
                            for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
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
                        results.add(new Tuple2<>((GeolocatedTS) ((Entry) c).getEntry(), intervals));
                }
            }

        } else {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance > epsilonSP)
                    continue;

                boolean found = false;
                for (int j = 0; j < GlobalConf.numClusters; j++) {

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);

                        int partition = Functions.getPartition(checkpoint);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;

                        int bundleIndex = (int) (checkpoint%(Math.ceil((double)GlobalConf.tsLength/GlobalConf.numTimePartitions)));
                        if ((queryTs[checkpoint] > c.getPartMaxTS()[partition][j][bundleIndex] && queryTs[checkpoint] - c.getPartMaxTS()[partition][j][bundleIndex] <= epsilonTS) ||
                                (queryTs[checkpoint] < c.getPartMinTS()[partition][j][bundleIndex] && c.getPartMinTS()[partition][j][bundleIndex] - queryTs[checkpoint] <= epsilonTS) ||
                                (queryTs[checkpoint] <= c.getPartMaxTS()[partition][j][bundleIndex] && queryTs[checkpoint] >= c.getPartMinTS()[partition][j][bundleIndex]) ) {

                            int partitionStart = GlobalConf.timePartitions[partition].get(0);
                            for (int k = checkpoint - 1; k >= partitionStart; k--) {
                                checks++;
                                bundleIndex = (int) (k%(Math.ceil((double)GlobalConf.tsLength/GlobalConf.numTimePartitions)));
                                if (queryTs[k] <= c.getPartMaxTS()[partition][j][bundleIndex] && queryTs[k] >= c.getPartMinTS()[partition][j][bundleIndex]) {
                                    start = k;
                                } else if (queryTs[k] > c.getPartMaxTS()[partition][j][bundleIndex]) {
                                    if (queryTs[k] - c.getPartMaxTS()[partition][j][bundleIndex] <= epsilonTS) {
                                        start = k;
                                    } else {
                                        break;
                                    }
                                } else if (queryTs[k] < c.getPartMinTS()[partition][j][bundleIndex]) {
                                    if (c.getPartMinTS()[partition][j][bundleIndex] - queryTs[k] <= epsilonTS) {
                                        start = k;
                                    } else {
                                        break;
                                    }
                                }
                            }

                            if (start == partitionStart) {
                                found = true;
                                break;
                            }

                            else {
                                partition = Functions.getPartition(checkpoint);
                                int partitionEnd = GlobalConf.timePartitions[partition].get(GlobalConf.timePartitions[partition].size() - 1);
                                for (int k = checkpoint + 1; k <= partitionEnd; k++) {
                                    checks++;
                                    bundleIndex = (int) (k % (Math.ceil((double) GlobalConf.tsLength / GlobalConf.numTimePartitions)));
                                    if (queryTs[k] <= c.getPartMaxTS()[partition][j][bundleIndex] && queryTs[k] >= c.getPartMinTS()[partition][j][bundleIndex]) {
                                        end = k;
                                    } else if (queryTs[k] > c.getPartMaxTS()[partition][j][bundleIndex]) {
                                        if (queryTs[k] - c.getPartMaxTS()[partition][j][bundleIndex] <= epsilonTS) {
                                            end = k;
                                        } else {
                                            break;
                                        }
                                    } else if (queryTs[k] < c.getPartMinTS()[partition][j][bundleIndex]) {
                                        if (c.getPartMinTS()[partition][j][bundleIndex] - queryTs[k] <= epsilonTS) {
                                            end = k;
                                        } else {
                                            break;
                                        }
                                    }
                                }

                                if ((end - start >= delta) || (end == partitionEnd)) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (found)
                        break;
                }
                if (found)
                    runQueryR(qcoords, queryTs, epsilonSP, epsilonTS, delta, c, results);

            }
        }
    }

    private static void runQueryR(double[] qcoords, double[] queryTs, double epsilonSP, double epsilonTS,
                                  int delta, Node n, ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) {
        nodeAccess++;
        if (n.leaf) {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance < epsilonSP) {
                    double[] ts = ((GeolocatedTS) ((Entry) c).getEntry()).getRawTimeSeries();

                    ArrayList<int[]> intervals = new ArrayList<>();

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;
                        if (Math.abs(ts[checkpoint] - queryTs[checkpoint]) <= epsilonTS) {
                            for (int j = checkpoint - 1; j >= 0; j--) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
                                    start = j;
                                } else {
                                    break;
                                }
                            }
                            for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) {
                                checks++;
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
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
                        results.add(new Tuple2<>((GeolocatedTS) ((Entry) c).getEntry(), intervals));
                }
            }

        } else {
            for (Node c : n.children) {
                double spatialDistance = Functions.distanceNode(c, qcoords);
                if (spatialDistance > epsilonSP)
                    continue;

                boolean found = false;
                for (int j = 0; j < c.getMaxTS().length; j++) {

                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint < end)
                            continue;
                        start = checkpoint;
                        end = checkpoint;

                        if ((queryTs[checkpoint] > c.getMaxTS()[j][checkpoint] && queryTs[checkpoint] - c.getMaxTS()[j][checkpoint] <= epsilonTS) ||
                                (queryTs[checkpoint] < c.getMinTS()[j][checkpoint] && c.getMinTS()[j][checkpoint] - queryTs[checkpoint] <= epsilonTS) ||
                                (queryTs[checkpoint] <= c.getMaxTS()[j][checkpoint] && queryTs[checkpoint] >= c.getMinTS()[j][checkpoint])) {

                            for (int k = checkpoint - 1; k >= 0; k--) {
                                checks++;
                                if (queryTs[k] <= c.getMaxTS()[j][k] && queryTs[k] >= c.getMinTS()[j][k]) {
                                    start = k;
                                } else if (queryTs[k] > c.getMaxTS()[j][k]) {
                                    if (queryTs[k] - c.getMaxTS()[j][k] <= epsilonTS) {
                                        start = k;
                                    } else {
                                        break;
                                    }
                                } else if (queryTs[k] < c.getMinTS()[j][k]) {
                                    if (c.getMinTS()[j][k] - queryTs[k] <= epsilonTS) {
                                        start = k;
                                    } else {
                                        break;
                                    }
                                }
                            }

                            for (int k = checkpoint + 1; k < GlobalConf.tsLength; k++) {
                                checks++;
                                if (queryTs[k] <= c.getMaxTS()[j][k] && queryTs[k] >= c.getMinTS()[j][k]) {
                                    end = k;
                                } else if (queryTs[k] > c.getMaxTS()[j][k]) {
                                    if (queryTs[k] - c.getMaxTS()[j][k] <= epsilonTS) {
                                        end = k;
                                    } else {
                                        break;
                                    }
                                } else if (queryTs[k] < c.getMinTS()[j][k]) {
                                    if (c.getMinTS()[j][k] - queryTs[k] <= epsilonTS) {
                                        end = k;
                                    } else {
                                        break;
                                    }
                                }
                            }

                            if (end - start >= delta) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (found)
                        break;
                }

                if (found)
                    runQueryR(qcoords, queryTs, epsilonSP, epsilonTS, delta, c, results);

            }
        }
    }
    */
}
