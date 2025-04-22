# First installation of spack 0.23.0

## Preparation

Start the following tutorial from a clean installation:

- Be careful with the content of your `.bashrc`, any environment variable or other package manager (conda, conan...) might get you into trouble
- Remove any old spack installation `$HOME/.spack`
- Start from a clean bash session: no modules loaded!

## Spack installation

### Download

```bash
wget https://github.com/spack/spack/releases/download/v0.23.0/spack-0.23.0.tar.gz
tar -xvf spack-0.23.0.tar.gz
rm spack-0.23.0.tar.gz
```

### Activation

```bash
. spack-0.23.0/share/spack/setup-env.sh
```

## Configuration

### Configure `config.yaml`

- Reduce the naming scheme of packages to avoid shebang issues.
- Increase the time out that is by default too short for some packages (like PDI).

```bash
spack config --scope site add 'config:connect_timeout:60'
```

### Configure `packages.yaml`

- Restrict `blas` and `lapack` providers to avoid `amdlibflame` and `amdblis`, for some reason they might get selected and fail to compile.

```bash
spack config --scope site add 'packages:all:providers:blas:[openblas]'
spack config --scope site add 'packages:all:providers:lapack:[openblas]'
```

### Configure `repos.yaml`

- Add PDI repository

```bash
git clone https://github.com/pdidev/spack.git spack-0.23.0/var/spack/repos/pdi
spack repo add --scope site spack-0.23.0/var/spack/repos/pdi
```

## Gyselalibxx environment installation

> :warning: This step takes some time

Using the following `gyselalibxx-env-0.23.0.yaml`

```yaml
spack:
  concretizer:
    unify: true
  definitions:
  - compilers:
    - 'gcc@11:'
  - packages:
    - 'benchmark'
    - 'cmake@3.22:3'
    - 'eigen'
    - 'fftw precision=float,double'
    - 'gdb ~debuginfod'
    - 'ginkgo@1.8'
    - 'googletest +gmock'
    - 'openmpi'
    - 'paraconf'
    - 'pdi@1.6:1'
    - 'pdiplugin-decl-hdf5 +mpi'
    - 'pdiplugin-decl-netcdf +mpi'
    - 'pdiplugin-mpi'
    - 'pdiplugin-set-value'
    - 'pdiplugin-trace'
  - python-packages:
    - 'py-dask'
    - 'py-h5netcdf'
    - 'py-h5py'
    - 'py-imageio'
    - 'py-matplotlib'
    - 'py-netcdf4 +mpi'
    - 'py-numpy'
    - 'py-scipy'
    - 'py-sympy'
    - 'py-xarray'
  specs:
  - lapack
  - matrix:
    - [$packages]
    - [$%compilers]
  - matrix:
    - [$python-packages]
    - [$%compilers]
  view: true
  packages:
    all:
      variants: ~mpi ~openmp cxxstd=17
```

```bash
spack env create gyselalibxx-env gyselalibxx-env-0.23.0.yaml
spack --env gyselalibxx-env install --jobs 2
```

## Enjoy

Ensure you start from a clean session

```bash
spack env activate -p gyselalibxx-env
```
