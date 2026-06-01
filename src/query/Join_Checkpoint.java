package queries;

import pbtsrindex.PBTSRTree;
import pbtsrindex.Entry;
import pbtsrindex.Node;
import tools.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;

//这是局部优化算法中的checkpoint方法

public class Join_Checkpoint {

    static ArrayList<Integer> checkpoints = new ArrayList<>();
    private static int nodeAccess;
    private static int checks;

    public Join_Checkpoint(){}

    //生成检查点
    private static void generateCheckpoints(int delta) {
        //delta是局部相似阈值
        int checkpoint = delta - 1;
        checkpoints.add(checkpoint);

        while (checkpoint < GlobalConf.tsLength) {
            checkpoint += delta;
            checkpoints.add(checkpoint);
        }
        checkpoints.remove(checkpoints.size() - 1);
    }

    //这里是设定不同查询方式和模式
    public static ArrayList<Tuple3<GeolocatedTS, GeolocatedTS, ArrayList<int[]>>> runJoin(PBTSRTree index1, PBTSRTree index2, double epsilonSP, double epsilonTS, int delta) throws IOException {
        LinkedList<Node> candidates1 = new LinkedList<>();
        LinkedList<Node> candidates2 = new LinkedList<>();
        ArrayList< Tuple3<GeolocatedTS, GeolocatedTS, ArrayList<int[]> > > results = new ArrayList<>();
        generateCheckpoints(delta); //生成checkpoint点

        for (Node child : index1.getRoot().children) {
            candidates1.add(child);
        }

        for (Node child : index2.getRoot().children) {
            candidates2.add(child);
        }

        runJoin(candidates1, candidates2, (index1 == index2), epsilonSP, epsilonTS, delta, results);
        GlobalConf.totalNodeAccess += nodeAccess;
        GlobalConf.totalChecks += checks;
        return results;
    }
   //这个另一个runjoin过程
    public static void runJoin(LinkedList<Node> candidates1, LinkedList<Node> candidates2, boolean isSelf, double epsilonSP, double epsilonTS, int delta,
                               ArrayList<Tuple3<GeolocatedTS, GeolocatedTS, ArrayList<int[]>>> results) throws IOException {

        for (int i = 0; i < candidates1.size(); i++) {
            for (int j = isSelf ? i : 0; j < candidates2.size(); j++) {

                LinkedList<Node> newCandidates1 = new LinkedList<>();
                LinkedList<Node> newCandidates2 = new LinkedList<>();

                //这里使用了bitvector算法
                if (!(candidates1.get(i).leaf) && !(candidates2.get(j).leaf)) {
                    Tuple2<Integer, Boolean> res = Functions.isInRangeNodeCheckpoint(candidates1.get(i), candidates2.get(j), epsilonSP, epsilonTS, delta, checkpoints, candidates1.get(i).getPartBitVector(),
                            candidates2.get(i).getPartBitVector());
                    checks += res.x;
                    if (res.y) {
                        nodeAccess += 2;
                        for (Node child : candidates1.get(i).children) {
                            newCandidates1.add(child);
                        }
                        for (Node child : candidates2.get(j).children) {
                            newCandidates2.add(child);
                        }
                        runJoin(newCandidates1, newCandidates2, candidates1.get(i) == candidates2.get(j), epsilonSP, epsilonTS, delta, results);
                    }
                } else if (candidates1.get(i).leaf && !(candidates1.get(i) instanceof Entry) && !(candidates2.get(j).leaf)) {
                    Tuple2<Integer, Boolean> res = Functions.isInRangeNodeCheckpoint(candidates1.get(i), candidates2.get(j), epsilonSP, epsilonTS, delta, checkpoints, candidates1.get(i).getPartBitVector(),
                            candidates2.get(i).getPartBitVector());
                    checks += res.x;
                    if (res.y) {
                        nodeAccess++;
                        newCandidates1.add(candidates1.get(i));
                        for (Node child : candidates2.get(j).children) {
                            newCandidates2.add(child);
                        }
                        runJoin(newCandidates1, newCandidates2, candidates1.get(i) == candidates2.get(j),
                                epsilonSP, epsilonTS, delta, results);
                    }
                } else if (!(candidates1.get(i).leaf) && !(candidates2.get(j) instanceof Entry) && candidates2.get(j).leaf) {
                    Tuple2<Integer, Boolean> res = Functions.isInRangeNodeCheckpoint(candidates1.get(i), candidates2.get(j), epsilonSP, epsilonTS, delta, checkpoints, candidates1.get(i).getPartBitVector(),
                            candidates2.get(i).getPartBitVector());
                    checks += res.x;
                    if (res.y) {
                        nodeAccess++;
                        for (Node child : candidates1.get(i).children) {
                            newCandidates1.add(child);
                        }
                        newCandidates2.add(candidates2.get(j));
                        runJoin(newCandidates1, newCandidates2, candidates1.get(i) == candidates2.get(j),
                                epsilonSP, epsilonTS, delta, results);
                    }
                } else if (candidates1.get(i).leaf && !(candidates1.get(i) instanceof Entry) && candidates2.get(j).leaf && !(candidates2.get(j) instanceof Entry)) {
                    Tuple2<Integer, Boolean> res = Functions.isInRangeNodeCheckpoint(candidates1.get(i), candidates2.get(j), epsilonSP, epsilonTS, delta, checkpoints, candidates1.get(i).getPartBitVector(),
                            candidates2.get(i).getPartBitVector());
                    checks += res.x;
                    if (res.y) {
                        nodeAccess += 2;

                        long startExecuteTime = System.currentTimeMillis();
                        GlobalConf.checks2++;
                        ArrayList<GeolocatedTS> contents1 = candidates1.get(i).getData();
                        ArrayList<GeolocatedTS> contents2 = candidates2.get(j).getData();
                        long totalElapsedExecuteTime = System.currentTimeMillis();
                        GlobalConf.tmpCount += totalElapsedExecuteTime - startExecuteTime;

                        if (candidates1.get(i) == candidates2.get(j)) {
                            for (int k = 0; k < contents1.size() - 1; k++) {
                                for (int n = k + 1; n < contents1.size(); n++) {
                                    double spatialDistance = Functions.euclideanDistance(contents1.get(k).getCoords(), contents1.get(n).getCoords());
                                    if (spatialDistance < epsilonSP) {
                                        double[] ts1 = contents1.get(k).getRawTimeSeries();
                                        double[] ts2 = contents1.get(n).getRawTimeSeries();
                                        ArrayList<int[]> intervals = checkTwoTS(ts1, ts2, epsilonTS, delta);
                                        if (intervals.size() > 0)
                                            results.add(new Tuple3<>(contents1.get(k), contents1.get(n), intervals));
                                    }
                                }
                            }
                        } else {
                            for (int k = 0; k < contents1.size(); k++) {
                                for (int n = 0; n < contents2.size(); n++) {
                                    double spatialDistance = Functions.euclideanDistance(contents1.get(k).getCoords(), contents2.get(n).getCoords());
                                    if (spatialDistance < epsilonSP) {
                                        double[] ts1 = contents1.get(k).getRawTimeSeries();
                                        double[] ts2 = contents2.get(n).getRawTimeSeries();
                                        ArrayList<int[]> intervals = checkTwoTS(ts1, ts2, epsilonTS, delta);
                                        if (intervals.size() > 0)
                                            results.add(new Tuple3<>(contents1.get(k), contents2.get(n), intervals));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
   // 论文中对sweep line 的优化算法  即实现了checkpoint方法
    private static ArrayList<int[]> checkTwoTS(double[] ts1, double[] ts2, double epsilonTS, int delta) {
        ArrayList<int[]> intervals = new ArrayList<>();
        int start, end = 0;
        for (int i = 0; i < checkpoints.size(); i++) {
            checks++;
            int checkpoint = checkpoints.get(i);
            if (checkpoint < end)
                continue;
            start = checkpoint;
            end = checkpoint;

            if (Math.abs(ts1[checkpoint] - ts2[checkpoint]) <= epsilonTS) {
                for (int j = checkpoint - 1; j >= 0; j--) {
                    checks++;
                    if (Math.abs(ts1[j] - ts2[j]) <= epsilonTS) {
                        start = j;
                    } else {
                        break;
                    }
                }
                for (int j = checkpoint + 1; j < GlobalConf.tsLength; j++) {
                    checks++;
                    if (Math.abs(ts1[j] - ts2[j]) <= epsilonTS) {
                        end = j;
                    } else {
                        break;
                    }
                }

                if (end - start >= delta)
                    intervals.add(new int[]{start, end});
            }
        }


        return intervals;
    }

}
