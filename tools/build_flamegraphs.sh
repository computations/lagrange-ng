for dir in */*/; do
  echo $dir
  pushd $dir
    if [ -f perf.data ]; then
      out_tmp=$(mktemp)
      perf script > $out_tmp
      stackcollapse-perf $out_tmp > out.folded.perf
      flamegraph out.folded.perf > out.svg
      rm $out_tmp
    fi
  popd
done
