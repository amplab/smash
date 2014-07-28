package edu.berkeley.cs.amplab.calldiff;

import java.util.Objects;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * A class encapsulating all the logic for {@link Object#hashCode} and {@link Object#equals}
 * methods. The client class should create a constant instance of this class and delegate the
 * {@code hashCode()} and {@code equals()} methods to it.
 */
public class HashCodeAndEquals<X> {

  @SafeVarargs
  public static <X> HashCodeAndEquals<X>
      create(Class<X> type, Function<? super X, ?>... accessors) {
    return new HashCodeAndEquals<>(
        obj -> Objects.hash(Stream.of(accessors)
            .map(
                new Function<Function<? super X, ?>, Object>() {
                  @Override public Object apply(Function<? super X, ?> accessor) {
                    return accessor.apply(obj);
                  }
                })
            .toArray()),
        (lhs, obj) -> {
          boolean same = lhs == obj;
          if (!same && type.isInstance(obj)) {
            X rhs = type.cast(obj);
            return Stream.of(accessors)
                .map(accessor -> Objects.equals(accessor.apply(lhs), accessor.apply(rhs)))
                .reduce(true, (x, y) -> x && y, (x, y) -> x && y);
          }
          return same;
        });
  }

  private final Function<X, Integer> hashCode;
  private final BiPredicate<X, Object> equals;

  private HashCodeAndEquals(Function<X, Integer> hashCode, BiPredicate<X, Object> equals) {
    this.hashCode = hashCode;
    this.equals = equals;
  }

  public int hashCode(X obj) {
    return hashCode.apply(obj);
  }

  public boolean equals(X lhs, Object obj) {
    return equals.test(lhs, obj);
  }
}
