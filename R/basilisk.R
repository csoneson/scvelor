#' @importFrom basilisk BasiliskEnvironment
scvelo <- basilisk::BasiliskEnvironment("scvelo",
                                        pkgname="scvelor",
                                        packages=c("numpy==1.18.4",
                                                   "scipy==1.4.1",
                                                   "numba==0.49.1",
                                                   "matplotlib==3.2.1",
                                                   "scikit-learn==0.23.0",
                                                   "h5py==2.10.0",
                                                   "click==7.1.2"),
                                        pip="scvelo==0.2.1")
