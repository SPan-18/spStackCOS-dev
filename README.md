# spStackCOS-dev

**An R package for Bayesian hierarchical regression of spatially-temporally misaligned outcome and exposure data.**

Spatial and temporal misalignment refers to the setting in which different variables are observed over incompatible spatial supports and/or at asynchronous time points or intervals. Such misalignment is common in studies of associations between human health outcomes and environmental exposures, such as air pollution indicators, which are often measured at different spatial and temporal resolutions. In this work, we develop

* A modular Bayesian inference framework that regresses an outcome on a spatially-temporally misaligned covariate
* Predictive stacking for analyzing spatially-temporally misaligned data

## ✨ Features

| Functionalities                                                                               | Supported |
|-----------------------------------------------------------------------------------------------|:---------:|
| Fit spatially point-referenced, temporally aggregated exposure data                           | ✅        |
| Perform predictive stacking for Bayesian inference                                            | ✅        |
| _Downscale_ posterior predictive inference at fine-scale space-time coordinates               | ✅        |
| _Upscale_ posterior predictive inference at polygons                                          | ✅        |
| Modular Bayesian inference for hierarchical regression                                        | ✅        |

## 🔧 Installation

```r
devtools::install_github("SPan-18/spStackCOS")
```