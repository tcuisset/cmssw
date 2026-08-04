# Mustache superclustering and its conditions

## Use in TICL

`TracksterLinkingbySuperClusteringMustache` considers tracksters in decreasing
raw transverse momentum. For each seed, a candidate must pass both:

1. `MustacheKernel::inDynamicDPhiWindow`, an energy-dependent maximum
   separation in phi; and
2. `MustacheKernel::inMustache`, the region between two parabolas in
   eta-phi space.

The calls use the seed trackster's eta and phi, but the candidate trackster's
raw energy, eta, and phi. Seed transverse energy is not an input to either
kernel.

## Dynamic delta-phi window

The `EcalSCDynamicDPhiParameters` payload contains parameter sets binned in
candidate energy and absolute seed eta. For a candidate,

```text
logEt = log10(candidateEnergy / cosh(candidateEta))
maxDPhi = yoffset + scale / (1 + exp((logEt - xoffset) / width))
maxDPhi = min(maxDPhi, cutoff)
maxDPhi = max(maxDPhi, saturation)
```

The candidate passes when the absolute, phi-wrapped separation from the seed
is smaller than `maxDPhi`.

## Mustache parabolas

`EcalMustacheSCParameters` contains:

- `sqrtLogClustETuning`;
- curvature coefficients `pUp` and `pLow`; and
- width coefficients `w0Up`, `w1Up`, `w0Low`, and `w1Low`.

The parameter set is selected using `log10(candidateEnergy)` and absolute
candidate eta. The coefficients produce an upper and lower boundary of the
form

```text
signedDeltaEta < upperCoefficient * deltaPhi^2 + upperOffset
signedDeltaEta > lowerCoefficient * deltaPhi^2 + lowerOffset
```

where `signedDeltaEta` is flipped for a negative-eta seed. The exact
implementation also applies a half-crystal-width term of `0.0087`. See
`RecoEcal/EgammaCoreTools/src/Mustache.cc` for the coefficient construction;
using `MustacheKernel` directly avoids duplicating its clipping and sign
conventions.

## Source of the parameters

In standard CMSSW reconstruction, the two ESProducts are conditions obtained
through the configured GlobalTag:

- `EcalMustacheSCParametersRcd`;
- `EcalSCDynamicDPhiParametersRcd`.

The GlobalTag maps each record to a payload tag, and the job's run selects the
applicable IOV. Consequently, the values must be quoted together with the
resolved GlobalTag and run. Run-2, Run-3, and Phase-2 GlobalTags can select
different optimized payloads.

CMSSW also contains
`EcalMustacheSCParametersESProducer_cff` and
`EcalSCDynamicDPhiParametersESProducer_cff`. They construct the objects from
hardcoded Python parameters for debugging, local replacement, and
CondTools/Ecal payload-writing workflows. Their presence does not mean they
are loaded by a standard `cmsDriver` configuration.

For example, `auto:phase2_realistic_T35` resolves in this release to
`150X_mcRun4_realistic_v1`, which maps the records to
`EcalMustacheSCParameters_average_mc` and
`EcalSCDynamicDPhiParameters_local_mc`. Those payload values happen to match
the parameter-producing Python fragments, but that is not guaranteed for
another GlobalTag.

## Dump and plot the active payload

`dumpMustacheParameters_cfg.py` loads the GlobalTag and schedules
`MustacheESProductDumper`. The analyzer consumes and prints the active typed
ESProducts, then calls the production kernels over an eta-phi grid:

```bash
cmsRun RecoHGCal/TICL/test/dumpMustacheParameters_cfg.py \
  seedEta=2.0 seedPhi=0.0 clusterEt=10 \
  output=mustache_scan.csv > mustache_parameters.log 2>&1

python3 RecoHGCal/TICL/test/plotMustache.py \
  mustache_scan.csv -o mustache_eta_phi.png
```

Use `clusterEt=-1 clusterEnergy=<raw-energy>` to hold candidate raw energy
fixed exactly as it is passed by the TICL linker.
