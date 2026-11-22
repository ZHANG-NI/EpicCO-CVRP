# EpicCO: Capacitated Vehicle Routing

The solver that [EpicCO](https://github.com/ZHANG-NI/EpicCO) produced for the capacitated vehicle routing problem, starting from HGS-CVRP (T. Vidal, Computers & Operations Research, 2022; https://github.com/vidalt/HGS-CVRP).

## Run

```
mkdir -p bin
bash build.sh bin
bin/solver <instance> <solution> -t <seconds> -seed <seed>
```

The runs of the paper add `-round 1` to the run command (rounded distances).
