package pl.genebeam.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;

public class RangeMultimapGeneral<R extends Comparable<R>, K> {
    private TreeRangeMap<R, List<K>> map = TreeRangeMap.create();

    public static void main(String[] args) {

        RangeMultimapGeneral<Integer, String> rangeMultimap = new RangeMultimapGeneral<Integer, String>();

        rangeMultimap.put(Range.closedOpen(0, 45), "0-45");
        rangeMultimap.put(Range.closedOpen(1, 5), "1-5");
        rangeMultimap.put(Range.closedOpen(3, 15), "3-15");
        rangeMultimap.put(Range.closedOpen(15, 20), "15-20");

        for (int i = -10; i < 50; i++) {
            System.out.println(i + "\t" + rangeMultimap.get(i));
        }

    }

    public Set<K> get(R key) {

        List<K> v = map.get(key);
        if (v == null)
            return Collections.emptySet();

        return new HashSet<K>(v);
    }

    public Set<K> get(Range<R> range) {
        RangeMap<R, List<K>> submap = map.subRangeMap(range);
        return submap.asMapOfRanges().values().stream().flatMap(Collection::stream).collect(Collectors.toSet());
    }

    public void put(Range<R> range, K value) {
        RangeMap<R, List<K>> subRangeMap = map.subRangeMap(range);
        Map<Range<R>, List<K>> submap = subRangeMap.asMapOfRanges();
        if (submap.isEmpty()) {
            ArrayList<K> list = new ArrayList<K>();
            list.add(value);
            map.put(range, list);
        } else {

            RangeSet<R> notCoveredSpan = TreeRangeSet.create();
            notCoveredSpan.add(range);

            List<Pair<Range<R>, List<K>>> newElements = new ArrayList<Pair<Range<R>, List<K>>>();
            for (Entry<Range<R>, List<K>> mapEntry : submap.entrySet()) {
                List<K> newValue = new ArrayList<K>(mapEntry.getValue());
                newValue.add(value);

                newElements.add(new ImmutablePair<Range<R>, List<K>>(mapEntry.getKey(), newValue));
                notCoveredSpan.remove(mapEntry.getKey());
            }

            for (Pair<Range<R>, List<K>> el : newElements) {
                map.put(el.getLeft(), el.getRight());
            }

            if (!notCoveredSpan.isEmpty()) {
                for (Range<R> notYetCoveredSpan : notCoveredSpan.asRanges()) {
                    ArrayList<K> list = new ArrayList<K>();
                    list.add(value);
                    map.put(notYetCoveredSpan, list);
                }
            }

        }
    }

}
