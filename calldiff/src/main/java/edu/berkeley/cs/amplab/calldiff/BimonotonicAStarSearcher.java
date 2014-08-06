package edu.berkeley.cs.amplab.calldiff;

import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.NONNULL;

import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.google.common.collect.AbstractIterator;

import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.RandomAccess;
import java.util.Set;
import java.util.Spliterators;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * A {@code BimonotonicAStarSearcher} is parameterized with 3 types: {@code X}, {@code Y}, and
 * {@code Z}, and is initialized with a {@link BiFunction} from {@code X} and {@code Y} to {@code Z}
 * and a {@link Comparator} on {@code Z}. The {@link #search} method accepts two
 * {@link RandomAccess} {@link List}s, one containing elements of type {@code X}, the other
 * containing elements of type {@code Y}. Each list is assumed to be in monotonically nondecreasing
 * order. The {@code BimonotonicAStarSearcher} will perform an A* search over the Cartesian product
 * of the two lists, using the {@code BiFunction} and {@code Comparator} to determine which element
 * to return next.
 */
public class BimonotonicAStarSearcher<X, Y, Z> {

  private static class AStar<X, Y> {

    static class Builder<X, Y> {

      private Comparator<? super X> comparator;
      private Function<? super X, ? extends Y> id;
      private Function<? super X, ? extends Stream<? extends X>> neighbors;

      AStar<X, Y> build() {
        return new AStar<>(comparator, id, neighbors);
      }

      Builder<X, Y> setComparator(Comparator<? super X> comparator) {
        this.comparator = comparator;
        return this;
      }

      Builder<X, Y> setIdFunction(Function<? super X, ? extends Y> id) {
        this.id = id;
        return this;
      }

      Builder<X, Y> setNeighborsFunction(
          Function<? super X, ? extends Stream<? extends X>> neighbors) {
        this.neighbors = neighbors;
        return this;
      }
    }

    static <X, Y> Builder<X, Y> builder() {
      return new Builder<>();
    }

    private final Comparator<? super X> comparator;
    private final Function<? super X, ? extends Y> id;
    private final Function<? super X, ? extends Stream<? extends X>> neighbors;

    private AStar(
        Comparator<? super X> comparator,
        Function<? super X, ? extends Y> id,
        Function<? super X, ? extends Stream<? extends X>> neighbors) {
      this.comparator = comparator;
      this.id = id;
      this.neighbors = neighbors;
    }

    Stream<X> search(X initial) {
      return StreamSupport.stream(
          Spliterators.spliteratorUnknownSize(
              new AbstractIterator<X>() {

                private final Queue<X> queue = new PriorityQueue<>(1, comparator);
                private final Set<Y> set = new HashSet<>();

                {
                  enqueue(initial);
                }

                @Override protected X computeNext() {
                  X next = queue.poll();
                  if (null != next) {
                    neighbors.apply(next).forEach(this::enqueue);
                    return next;
                  }
                  return endOfData();
                }

                private void enqueue(X x) {
                  if (set.add(id.apply(x))) {
                    queue.offer(x);
                  }
                }
              },
              IMMUTABLE | NONNULL),
          false);
    }
  }

  public static class Builder<X, Y, Z> {

    private BiFunction<? super X, ? super Y, ? extends Z> biFunction;
    private Comparator<? super Z> comparator;

    private Builder() {}

    public BimonotonicAStarSearcher<X, Y, Z> build() {
      return new BimonotonicAStarSearcher<>(biFunction, comparator);
    }

    public Builder<X, Y, Z> setBiFunction(
        BiFunction<? super X, ? super Y, ? extends Z> biFunction) {
      this.biFunction = biFunction;
      return this;
    }

    public Builder<X, Y, Z> setComparator(Comparator<? super Z> comparator) {
      this.comparator = comparator;
      return this;
    }
  }

  private static class Pair {

    private static final HashCodeAndEquals<Pair> HASH_CODE_AND_EQUALS = HashCodeAndEquals.create(
        Pair.class,
        Pair::i,
        Pair::j);

    private final int i, j;

    Pair(int i, int j) {
      this.i = i;
      this.j = j;
    }

    @Override public boolean equals(Object obj) {
      return HASH_CODE_AND_EQUALS.equals(this, obj);
    }

    @Override public int hashCode() {
      return HASH_CODE_AND_EQUALS.hashCode(this);
    }

    int i() {
      return i;
    }

    int j() {
      return j;
    }
  }

  public static <X, Y, Z> Builder<X, Y, Z> builder() {
    return new Builder<>();
  }

  private final BiFunction<? super X, ? super Y, ? extends Z> biFunction;
  private final Comparator<? super Z> comparator;

  private BimonotonicAStarSearcher(
      BiFunction<? super X, ? super Y, ? extends Z> biFunction,
      Comparator<? super Z> comparator) {
    this.biFunction = biFunction;
    this.comparator = comparator;
  }

  public <L extends List<? extends X> & RandomAccess, R extends List<? extends Y> & RandomAccess>
      Stream<Z> search(L lhs, R rhs) {
    int lhsSize = lhs.size(), rhsSize = rhs.size();
    class PairWithValue {

      private final Pair pair;
      private final Supplier<Z> value;

      PairWithValue(int i, int j) {
        this.pair = new Pair(i, j);
        this.value = Suppliers.memoize(
            () -> biFunction.apply(lhs.get(pair.i()), rhs.get(pair.j())));
      }

      Stream<PairWithValue> neighbors() {
        int i = pair.i(), nextI = i + 1, j = pair.j(), nextJ = j + 1;
        return Stream.concat(
            nextI < lhsSize ? Stream.of(new PairWithValue(nextI, j)) : Stream.empty(),
            nextJ < rhsSize ? Stream.of(new PairWithValue(i, nextJ)) : Stream.empty());
      }

      Pair pair() {
        return pair;
      }

      Z value() {
        return value.get();
      }
    }
    return 0 == lhsSize || 0 == rhsSize
        ? Stream.empty()
        : AStar.<PairWithValue, Pair>builder()
            .setComparator(Comparator.comparing(PairWithValue::value, this.comparator))
            .setIdFunction(PairWithValue::pair)
            .setNeighborsFunction(PairWithValue::neighbors)
            .build()
            .search(new PairWithValue(0, 0))
            .map(PairWithValue::value);
  }
}
