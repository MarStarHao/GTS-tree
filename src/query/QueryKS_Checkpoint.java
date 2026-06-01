package queries;

import com.google.common.collect.MinMaxPriorityQueue;
import pbtsrindex.Entry;
import pbtsrindex.Node;
import tools.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.PriorityQueue;
/*
* 主要是
* */

public class QueryKS_Checkpoint {

    static ArrayList<Integer> checkpoints;
    private static int nodeAccess;
    private static int checks;

    private static void generateCheckpoints(int delta) {
        checkpoints = new ArrayList<>();
        int checkpoint = delta - 1;
        checkpoints.add(checkpoint);

        while (checkpoint < GlobalConf.tsLength) {
            checkpoint += delta;
            checkpoints.add(checkpoint);
        }
        checkpoints.remove(checkpoints.size() - 1);
    }

    public static ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> runQueryKS(double[] coords,
                                                                               double[] queryTs,
                                                                               double epsilonSP,
                                                                               double epsilonTS,
                                                                               int k, int delta,
                                                                               Node root) throws IOException {

        generateCheckpoints(delta); //生成一个检查点
        ArrayList< Tuple2< GeolocatedTS, ArrayList<int[]> > > results = new ArrayList<>(); //得到一个二元组
        nodeAccess = 0;
        runQueryKS(coords, queryTs, epsilonSP, epsilonTS, k, delta, root, results);
        GlobalConf.totalNodeAccess += nodeAccess;
        GlobalConf.totalChecks += checks;
        return results;
    }

    public static ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> runQueryKS(double[] qcoords, double[] queryTs, double epsilonSP, double epsilonTS, int kappa,
                                                                               int delta, Node node, ArrayList<Tuple2<GeolocatedTS, ArrayList<int[]>>> results) throws IOException {

        PriorityQueue<IntervalSize> queueNode = new PriorityQueue<>();  //kappa是参数k
        MinMaxPriorityQueue<IntervalSize> queueRaw = MinMaxPriorityQueue.orderedBy(new DistComparator()).maximumSize(kappa).create();
        queueNode.add(new IntervalSize(node, new int[]{0, 0}, 0));
        nodeAccess = 0;

        while (queueNode.isEmpty() == false) {
            IntervalSize ed = queueNode.poll();
            if (queueRaw.size() > 0) {
                if (ed.size < queueRaw.peekLast().size) {
                    break;
                }
            }
            Node n = ed.n;
            if (n.leaf && !(n instanceof Entry)) {
                nodeAccess++;
                ArrayList<GeolocatedTS> contents = n.getData();
                for (GeolocatedTS geoTS : contents) {
                    double spatialDistance = Functions.euclideanDistance(geoTS.getCoords(), qcoords);
                    if (spatialDistance > epsilonSP)
                        continue;

                    double[] ts = geoTS.getRawTimeSeries();
                    int[] largestInterval = new int[]{0, 0};
                    int start, end = 0;
                    for (int i = 0; i < checkpoints.size(); i++) {
                        checks++;

                        int checkpoint = checkpoints.get(i);
                        if (checkpoint <= end)
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
                                if (Math.abs(ts[j] - queryTs[j]) <= epsilonTS) {
                                    end = j;
                                    if (j == GlobalConf.tsLength - 1) {
                                        if (end - start >= delta) {
                                            if (end - start > largestInterval[1] - largestInterval[0])
                                                largestInterval = new int[]{start, end};
                                        }
                                    }
                                } else {
                                    if (end - start >= delta) {
                                        if (end - start > largestInterval[1] - largestInterval[0])
                                            largestInterval = new int[]{start, end};
                                    }
                                    break;
                                }
                            }

                            if (end - start >= delta) {
                                if (end - start > largestInterval[1] - largestInterval[0])
                                    largestInterval = new int[]{start, end};
                            }
                        }
                    }

                    int dur = largestInterval[1] - largestInterval[0];
                    if (queueRaw.size() < kappa) {
                        queueRaw.add(new IntervalSize(new Entry<>(geoTS.getCoords(), new double[]{1, 1}, geoTS), largestInterval, dur));
                    } else {
                        queueRaw.add(new IntervalSize(new Entry<>(geoTS.getCoords(), new double[]{1, 1}, geoTS), largestInterval, dur));
                        delta = queueRaw.peekLast().size;
                    }
                }
                generateCheckpoints(delta);
            } else {
                nodeAccess++;
                for (Node c : n.children) {

                    double spatialDistance = Functions.distanceNode(c, qcoords);
                    if (spatialDistance > epsilonSP)
                        continue;

                    if (!c.leaf) {
                        queueNode.add(new IntervalSize(c, new int[]{0, GlobalConf.tsLength}, GlobalConf.tsLength));
                    } else {
                        int[][][] bitVector = c.getPartBitVector();
                        Tuple2<Integer, int[]> res = Functions.getLargestInterval(queryTs, c, checkpoints, epsilonTS, delta, bitVector);
                        int dur = res.y[1] - res.y[0];
                        if (queueRaw.size() < kappa) {
                            queueNode.add(new IntervalSize(c, res.y, dur));
                        }
                        else {
                            if (dur >= delta) {
                                queueNode.add(new IntervalSize(c, res.y, dur));
                            }
                        }
                    }
                }
            }
        }

        for (int i=0; i<kappa; i++) {
            IntervalSize el = queueRaw.pollFirst();
            ArrayList<int[]> intervals = new ArrayList<>();
            intervals.add(el.interval);
            results.add(new Tuple2<>((GeolocatedTS) ((Entry) el.n).getEntry(), intervals));
        }

        return results;
    }
}