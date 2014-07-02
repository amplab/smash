package edu.berkeley.cs.amplab.smash4j;

import com.google.common.base.Function;
import com.google.common.base.Functions;
import com.google.common.base.Optional;
import com.google.common.base.Preconditions;
import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;
import com.google.common.collect.FluentIterable;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Iterables;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;

import edu.berkeley.cs.amplab.smash4j.vcf.VcfRecord;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public abstract class Variant {

  public static class Builder {

    private List<String> alternateBases;
    private List<Call> calls;
    private String contig;
    private ListMultimap<String, String> info;
    private List<String> names;
    private Integer position;
    private String referenceBases;

    private Builder(
        List<String> alternateBases,
        List<Call> calls,
        String contig,
        ListMultimap<String, String> info,
        List<String> names,
        Integer position,
        String referenceBases) {
      if (null != alternateBases) {
        setAlternateBases(alternateBases);
      }
      if (null != calls) {
        setCalls(calls);
      }
      if (null != contig) {
        setContig(contig);
      }
      if (null != info) {
        setInfo(info);
      }
      if (null != names) {
        setNames(names);
      }
      if (null != position) {
        setPosition(position);
      }
      if (null != referenceBases) {
        setReferenceBases(referenceBases);
      }
    }

    public Variant build() {
      Preconditions.checkNotNull(alternateBases);
      Preconditions.checkNotNull(calls);
      Preconditions.checkNotNull(contig);
      Preconditions.checkNotNull(info);
      Preconditions.checkNotNull(names);
      Preconditions.checkNotNull(position);
      Preconditions.checkNotNull(referenceBases);
      return new Variant() {

            @Override public List<String> alternateBases() {
              return alternateBases;
            }

            @Override public List<Call> calls() {
              return calls;
            }

            @Override public String contig() {
              return contig;
            }

            @Override public ListMultimap<String, String> info() {
              return info;
            }

            @Override public List<String> names() {
              return names;
            }

            @Override public int position() {
              return position;
            }

            @Override public String referenceBases() {
              return referenceBases;
            }
          };
    }

    public Builder setAlternateBases(List<String> alternateBases) {
      this.alternateBases = alternateBases;
      return this;
    }

    public Builder setCalls(List<Call> calls) {
      this.calls = calls;
      return this;
    }

    public Builder setContig(String contig) {
      this.contig = contig;
      return this;
    }

    public Builder setInfo(ListMultimap<String, String> info) {
      this.info = info;
      return this;
    }

    public Builder setNames(List<String> names) {
      this.names = names;
      return this;
    }

    public Builder setPosition(int position) {
      this.position = position;
      return this;
    }

    public Builder setReferenceBases(String referenceBases) {
      this.referenceBases = referenceBases;
      return this;
    }
  }

  public abstract static class Call {

    public static class Builder {

      private List<Integer> genotype;
      private List<Double> genotypeLikelihood;
      private ListMultimap<String, String> info;
      private Optional<Integer> phaseset;

      private Builder(
          List<Integer> genotype,
          List<Double> genotypeLikelihood,
          ListMultimap<String, String> info,
          Optional<Integer> phaseset) {
        if (null != genotype) {
          setGenotype(genotype);
        }
        if (null != genotypeLikelihood) {
          setGenotypeLikelihood(genotypeLikelihood);
        }
        if (null != info) {
          setInfo(info);
        }
        if (null != phaseset) {
          setPhaseset(phaseset);
        }
      }

      public Call build() {
        Preconditions.checkNotNull(genotype);
        Preconditions.checkNotNull(genotypeLikelihood);
        Preconditions.checkNotNull(info);
        Preconditions.checkNotNull(phaseset);
        return new Call() {

              @Override public List<Integer> genotype() {
                return genotype;
              }

              @Override public List<Double> genotypeLikelihood() {
                return genotypeLikelihood;
              }

              @Override public ListMultimap<String, String> info() {
                return info;
              }

              @Override public Optional<Integer> phaseset() {
                return phaseset;
              }
            };
      }

      public Builder setGenotype(List<Integer> genotype) {
        this.genotype = genotype;
        return this;
      }

      public Builder setGenotypeLikelihood(List<Double> genotypeLikelihood) {
        this.genotypeLikelihood = genotypeLikelihood;
        return this;
      }

      public Builder setInfo(ListMultimap<String, String> info) {
        this.info = info;
        return this;
      }

      public Builder setPhaseset(Optional<Integer> phaseset) {
        this.phaseset = phaseset;
        return this;
      }
    }

    static final Function<com.google.api.services.genomics.model.Call, Call> CREATE_FROM_CALL =
        new Function<com.google.api.services.genomics.model.Call, Call>() {
          @Override public Call apply(com.google.api.services.genomics.model.Call call) {
            return create(call);
          }
        };

    public static final Function<Call, List<Integer>> GENOTYPE =
        new Function<Call, List<Integer>>() {
          @Override public List<Integer> apply(Call call) {
            return call.genotype();
          }
        };

    public static final Function<Call, List<Double>> GENOTYPE_LIKELIHOOD =
        new Function<Call, List<Double>>() {
          @Override public List<Double> apply(Call call) {
            return call.genotypeLikelihood();
          }
        };

    public static final Function<Call, ListMultimap<String, String>> INFO =
        new Function<Call, ListMultimap<String, String>>() {
          @Override public ListMultimap<String, String> apply(Call call) {
            return call.info();
          }
        };

    public static final Function<Call, Optional<Integer>> PHASESET =
        new Function<Call, Optional<Integer>>() {
          @Override public Optional<Integer> apply(Call call) {
            return call.phaseset();
          }
        };

    public static Builder builder() {
      return new Builder(
          Collections.<Integer>emptyList(),
          Collections.<Double>emptyList(),
          ImmutableListMultimap.<String, String>of(),
          Optional.<Integer>absent());
    }

    static Call create(final com.google.api.services.genomics.model.Call delegate) {
      return new Call() {

        private final Supplier<List<Integer>> genotype = Suppliers.memoize(
            new Supplier<List<Integer>>() {
              @Override public List<Integer> get() {
                return Optional.fromNullable(delegate.getGenotype())
                    .or(Collections.<Integer>emptyList());
              }
            });

        private final Supplier<List<Double>> genotypeLikelihood = Suppliers.memoize(
            new Supplier<List<Double>>() {
              @Override public List<Double> get() {
                return Optional.fromNullable(delegate.getGenotypeLikelihood())
                    .or(Collections.<Double>emptyList());
              }
            });

        private final Supplier<ListMultimap<String, String>> info = Suppliers.memoize(
            new Supplier<ListMultimap<String, String>>() {
              @Override public ListMultimap<String, String> get() {
                return toMultimap(Optional.fromNullable(delegate.getInfo())
                    .or(Collections.<String, List<String>>emptyMap()));
              }
            });

        private final Supplier<Optional<Integer>> phaseset = Suppliers.memoize(
            new Supplier<Optional<Integer>>() {
              @Override public Optional<Integer> get() {
                return Optional.fromNullable(delegate.getPhaseset()).transform(PARSE_INT);
              }
            });

        @Override public List<Integer> genotype() {
          return genotype.get();
        }

        @Override public List<Double> genotypeLikelihood() {
          return genotypeLikelihood.get();
        }

        @Override public ListMultimap<String, String> info() {
          return info.get();
        }

        @Override public Optional<Integer> phaseset() {
          return phaseset.get();
        }
      };
    }

    static Call create(final List<String> format, final List<String> sample) {
      return new Call() {

        private final Supplier<List<Integer>> genotype = Suppliers.memoize(
            new Supplier<List<Integer>>() {
              @Override public List<Integer> get() {
                return GENOTYPE_SPLITTER.apply(
                    Iterables.getOnlyElement(originalInfo.get().get("GT")));
              }
            });

        private final Supplier<List<Double>> genotypeLikelihood = Suppliers.memoize(
            new Supplier<List<Double>>() {
              @Override public List<Double> get() {
                return FluentIterable.from(originalInfo.get().get("GL"))
                    .transform(PARSE_DOUBLE)
                    .toList();
              }
            });

        private final Supplier<ListMultimap<String, String>> info = Suppliers.memoize(
            new Supplier<ListMultimap<String, String>>() {
              @Override public ListMultimap<String, String> get() {
                ImmutableListMultimap.Builder<String, String>
                    multimap = ImmutableListMultimap.builder();
                for (Map.Entry<String, String> entry : originalInfo.get().entries()) {
                  String key = entry.getKey();
                  switch (key) {
                    case "GT":
                    case "GL":
                    case "PS":
                      break;
                    default:
                      multimap.put(key, entry.getValue());
                  }
                }
                return multimap.build();
              }
            });

        private final Supplier<ListMultimap<String, String>> originalInfo = Suppliers.memoize(
            new Supplier<ListMultimap<String, String>>() {

              @Override public ListMultimap<String, String> get() {
                return toMultimap(Maps.transformValues(zip(format, sample), COMMA_SPLITTER));
              }

              private <K, V> Map<K, V>
                  zip(Iterable<? extends K> keys, Iterable<? extends V> values) {
                Iterator<? extends K> keysIterator = keys.iterator();
                Iterator<? extends V> valuesIterator = values.iterator();
                ImmutableMap.Builder<K, V> map = ImmutableMap.builder();
                for (; keysIterator.hasNext() && valuesIterator.hasNext();
                    map.put(keysIterator.next(), valuesIterator.next()));
                Preconditions.checkState(!(keysIterator.hasNext() || valuesIterator.hasNext()));
                return map.build();
              }
            });

        private final Supplier<Optional<Integer>> phaseset = Suppliers.memoize(
            new Supplier<Optional<Integer>>() {
              @Override public Optional<Integer> get() {
                List<String> phaseset = originalInfo.get().get("PS");
                switch (phaseset.size()) {
                  case 0:
                    return Optional.absent();
                  case 1:
                    return Optional.of(Iterables.getOnlyElement(phaseset)).transform(PARSE_INT);
                  default:
                    throw new IllegalStateException();
                }
              }
            });

        @Override public List<Integer> genotype() {
          return genotype.get();
        }

        @Override public List<Double> genotypeLikelihood() {
          return genotypeLikelihood.get();
        }

        @Override public ListMultimap<String, String> info() {
          return info.get();
        }

        @Override public Optional<Integer> phaseset() {
          return phaseset.get();
        }
      };
    }

    @Override
    public final boolean equals(Object obj) {
      boolean same = this == obj;
      if (!same && obj instanceof Call) {
        Call rhs = (Call) obj;
        return Objects.equals(genotype(), rhs.genotype())
            && Objects.equals(genotypeLikelihood(), rhs.genotypeLikelihood())
            && Objects.equals(info(), rhs.info())
            && Objects.equals(phaseset(), rhs.phaseset());
      }
      return same;
    }

    public abstract List<Integer> genotype();

    public abstract List<Double> genotypeLikelihood();

    @Override
    public final int hashCode() {
      return Objects.hash(
          genotype(),
          genotypeLikelihood(),
          info(),
          phaseset());
    }

    public abstract ListMultimap<String, String> info();

    public abstract Optional<Integer> phaseset();

    public final Builder toBuilder() {
      return new Builder(genotype(), genotypeLikelihood(), info(), phaseset());
    }

    @Override
    public final String toString() {
      return String.format(
          "genotype: %s "
              + "genotypeLikelihood: %s "
              + "info: %s "
              + "phaseset: %s",
          genotype(),
          genotypeLikelihood(),
          info(),
          phaseset());
    }
  }

  public static final Function<Variant, List<String>> ALTERNATE_BASES =
      new Function<Variant, List<String>>() {
        @Override public List<String> apply(Variant variant) {
          return variant.alternateBases();
        }
      };

  public static final Function<Variant, List<Call>> CALLS =
      new Function<Variant, List<Call>>() {
        @Override public List<Call> apply(Variant variant) {
          return variant.calls();
        }
      };

  private static final Function<CharSequence, List<String>> COMMA_SPLITTER =
      splitter("([^,]+?)(?:,|$)", Functions.<String>identity());

  public static final Function<Variant, String> CONTIG =
      new Function<Variant, String>() {
        @Override public String apply(Variant variant) {
          return variant.contig();
        }
      };

  public static final Function<com.google.api.services.genomics.model.Variant, Variant>
      CREATE_FROM_VARIANT =
      new Function<com.google.api.services.genomics.model.Variant, Variant>() {
        @Override public Variant apply(com.google.api.services.genomics.model.Variant variant) {
          return create(variant);
        }
      };

  public static final Function<VcfRecord, Variant> CREATE_FROM_VCF_RECORD =
      new Function<VcfRecord, Variant>() {
        @Override public Variant apply(VcfRecord record) {
          return create(record);
        }
      };

  private static final Function<String, Integer> PARSE_INT =
      new Function<String, Integer>() {
        @Override public Integer apply(String input) {
          return Integer.parseInt(input);
        }
      };

  static final Function<CharSequence, List<Integer>> GENOTYPE_SPLITTER =
      splitter("(0|(?:[1-9][0-9]*?))(?:[/|]|$)", PARSE_INT);

  public static final Function<Variant, ListMultimap<String, String>> INFO =
      new Function<Variant, ListMultimap<String, String>>() {
        @Override public ListMultimap<String, String> apply(Variant variant) {
          return variant.info();
        }
      };

  public static final Function<Variant, List<String>> NAMES =
      new Function<Variant, List<String>>() {
        @Override public List<String> apply(Variant variant) {
          return variant.names();
        }
      };

  private static final Function<String, Double> PARSE_DOUBLE =
      new Function<String, Double>() {
        @Override public Double apply(String input) {
          return Double.parseDouble(input);
        }
      };

  public static final Function<Variant, Integer> POSITION =
      new Function<Variant, Integer>() {
        @Override public Integer apply(Variant variant) {
          return variant.position();
        }
      };

  public static final Function<Variant, String> REFERENCE_BASES =
      new Function<Variant, String>() {
        @Override public String apply(Variant variant) {
          return variant.referenceBases();
        }
      };

  public static Builder builder() {
    return new Builder(
        Collections.<String>emptyList(),
        Collections.<Call>emptyList(),
        null,
        ImmutableListMultimap.<String, String>of(),
        Collections.<String>emptyList(),
        null,
        null);
  }

  public static Variant create(final com.google.api.services.genomics.model.Variant delegate) {
    return new Variant() {

      private final Supplier<List<String>> alternateBases = Suppliers.memoize(
          new Supplier<List<String>>() {
            @Override public List<String> get() {
              return Optional.fromNullable(delegate.getAlternateBases())
                  .or(Collections.<String>emptyList());
            }
          });

      private final Supplier<List<Call>> calls = Suppliers.memoize(
          new Supplier<List<Call>>() {
            @Override public List<Call> get() {
              return FluentIterable
                  .from(
                      Optional.fromNullable(delegate.getCalls())
                          .or(Collections.<com.google.api.services.genomics.model.Call>emptyList()))
                  .transform(Call.CREATE_FROM_CALL)
                  .toList();
            }
          });

      private final Supplier<String> contig = Suppliers.memoize(
          new Supplier<String>() {
            @Override public String get() {
              return delegate.getContig();
            }
          });

      private final Supplier<ListMultimap<String, String>> info = Suppliers.memoize(
          new Supplier<ListMultimap<String, String>>() {
            @Override public ListMultimap<String, String> get() {
              return toMultimap(
                  Optional.fromNullable(delegate.getInfo())
                      .or(Collections.<String, List<String>>emptyMap()));
            }
          });

      private final Supplier<List<String>> names = Suppliers.memoize(
          new Supplier<List<String>>() {
            @Override public List<String> get() {
              return Optional.fromNullable(delegate.getNames())
                  .or(Collections.<String>emptyList());
            }
          });

      private final Supplier<Integer> position = Suppliers.memoize(
          new Supplier<Integer>() {
            @Override public Integer get() {
              return delegate.getPosition().intValue();
            }
          });

      private final Supplier<String> referenceBases = Suppliers.memoize(
          new Supplier<String>() {
            @Override public String get() {
              return delegate.getReferenceBases();
            }
          });

      @Override public List<String> alternateBases() {
        return alternateBases.get();
      }

      @Override public List<Call> calls() {
        return calls.get();
      }

      @Override public String contig() {
        return contig.get();
      }

      @Override public ListMultimap<String, String> info() {
        return info.get();
      }

      @Override public List<String> names() {
        return names.get();
      }

      @Override public int position() {
        return position.get();
      }

      @Override public String referenceBases() {
        return referenceBases.get();
      }
    };
  }

  public static Variant create(final VcfRecord delegate) {
    return new Variant() {

      private final Supplier<List<String>> alternateBases = Suppliers.memoize(
          new Supplier<List<String>>() {
            @Override public List<String> get() {
              return Optional.fromNullable(delegate.alt())
                  .or(Collections.<String>emptyList());
            }
          });

      private final Supplier<List<Call>> calls = Suppliers.memoize(
          new Supplier<List<Call>>() {
            @Override
            public List<Call> get() {
              return FluentIterable
                  .from(
                      Optional.fromNullable(delegate.samples())
                          .or(Collections.<List<String>>emptyList()))
                  .transform(
                      new Function<List<String>, Call>() {

                        private final List<String> format =
                            Optional.fromNullable(delegate.format())
                                .or(Collections.<String>emptyList());

                        @Override public Call apply(final List<String> sample) {
                          return Call.create(format, sample);
                        }
                      })
                  .toList();
            }
          });

      private final Supplier<String> contig = Suppliers.memoize(
          new Supplier<String>() {
            @Override public String get() {
              return delegate.chrom();
            }
          });

      private final Supplier<ListMultimap<String, String>> info = Suppliers.memoize(
          new Supplier<ListMultimap<String, String>>() {
            @Override public ListMultimap<String, String> get() {
              return toMultimap(
                  Optional.fromNullable(delegate.info())
                      .or(Collections.<String, List<String>>emptyMap()));
            }
          });

      private final Supplier<List<String>> names = Suppliers.memoize(
          new Supplier<List<String>>() {
            @Override public List<String> get() {
              return Optional.fromNullable(delegate.ids())
                  .or(Collections.<String>emptyList());
            }
          });

      private final Supplier<Integer> position = Suppliers.memoize(
          new Supplier<Integer>() {
            @Override public Integer get() {
              return delegate.pos();
            }
          });

      private final Supplier<String> referenceBases = Suppliers.memoize(
          new Supplier<String>() {
            @Override public String get() {
              return delegate.ref();
            }
          });

      @Override public List<String> alternateBases() {
        return alternateBases.get();
      }

      @Override
      public List<Call> calls() {
        return calls.get();
      }

      @Override public String contig() {
        return contig.get();
      }

      @Override public ListMultimap<String, String> info() {
        return info.get();
      }

      @Override public List<String> names() {
        return names.get();
      }

      @Override public int position() {
        return position.get();
      }

      @Override public String referenceBases() {
        return referenceBases.get();
      }
    };
  }

  private static <X> Function<CharSequence, List<X>> splitter(
      final String regex, final Function<? super String, ? extends X> transformer) {
    return new Function<CharSequence, List<X>>() {

          private final Pattern pattern = Pattern.compile(regex);

          @Override public List<X> apply(CharSequence input) {
            ImmutableList.Builder<X> list = ImmutableList.builder();
            for (Matcher matcher = pattern.matcher(input); matcher.find();
                list.add(transformer.apply(matcher.group(1))));
            return list.build();
          }
        };
  }

  private static <X, Y> ListMultimap<X, Y> toMultimap(Map<X, List<Y>> map) {
    ImmutableListMultimap.Builder<X, Y> multimap = ImmutableListMultimap.builder();
    for (Map.Entry<X, List<Y>> entry : map.entrySet()) {
      multimap.putAll(entry.getKey(), entry.getValue());
    }
    return multimap.build();
  }

  public abstract List<String> alternateBases();

  public abstract List<Call> calls();

  public abstract String contig();

  @Override public final boolean equals(Object obj) {
    boolean same = this == obj;
    if (!same && obj instanceof Variant) {
      Variant rhs = (Variant) obj;
      return Objects.equals(alternateBases(), rhs.alternateBases())
          && Objects.equals(calls(), rhs.calls()) && Objects.equals(contig(), rhs.contig())
          && Objects.equals(info(), rhs.info()) && Objects.equals(names(), rhs.names())
          && Objects.equals(position(), rhs.position())
          && Objects.equals(referenceBases(), rhs.referenceBases());
    }
    return same;
  }

  @Override public final int hashCode() {
    return Objects.hash(alternateBases(),
        calls(),
        contig(),
        info(),
        names(),
        position(),
        referenceBases());
  }

  public abstract ListMultimap<String, String> info();

  public abstract List<String> names();

  public abstract int position();

  public abstract String referenceBases();

  public final Builder toBuilder() {
    return new Builder(
        alternateBases(),
        calls(),
        contig(),
        info(),
        names(),
        position(),
        referenceBases());
  }

  @Override public final String toString() {
    return String.format(
        "alternateBases: %s "
            + "calls: %s "
            + "contig: %s "
            + "info: %s "
            + "names: %s "
            + "position: %s "
            + "referenceBases: %s",
        alternateBases(),
        calls(),
        contig(),
        info(),
        names(),
        position(),
        referenceBases());
  }
}
