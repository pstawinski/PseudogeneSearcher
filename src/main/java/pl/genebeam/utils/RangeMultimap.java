package pl.genebeam.utils;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

public class RangeMultimap<K> extends RangeMultimapGeneral<Integer, K> {
    private TreeRangeMap<Integer, List<K>> map = TreeRangeMap.create();

    public static void main(String[] args) {

        RangeMultimap<String> rangeMultimap = new RangeMultimap<String>();

        rangeMultimap.put(Range.closedOpen(0, 45), "0-45");
        rangeMultimap.put(Range.closedOpen(1, 5), "1-5");
        rangeMultimap.put(Range.closedOpen(3, 15), "3-15");
        rangeMultimap.put(Range.closedOpen(15, 20), "15-20");

        for (int i = -10; i < 50; i++) {
            System.out.println(i + "\t" + rangeMultimap.get(i));
        }

        for (int i = -10; i < 50; i++) {
            System.out.println(i + "\t" + rangeMultimap.get(i, 3));
        }
    }

    public Set<K> get(int key, int radius) {
        if (radius == 1) {
            return get(key);
        }

        RangeMap<Integer, List<K>> submap = map.subRangeMap(Range.closed(key - radius, key + radius));

        Set<K> output = new HashSet<K>();
        for (List<K> values : submap.asMapOfRanges().values()) {
            output.addAll(values);
        }

        return output;
    }

}
