Name:           wannier90
Summary:        Maximally-Localised Generalised Wannier Functions Code
Version:        0.0.0
Release:        %autorelease
License:        GPLv2
URL:            https://www.wannier.org/

Source:         https://github.com/wannier-developers/wannier90/archive/refs/tags/v%{version}.tar.gz

BuildRequires:  ninja-build
BuildRequires:  cmake
BuildRequires:  gcc-fortran
BuildRequires:  flexiblas-devel
# Required for testing
BuildRequires:  gcc-c++
BuildRequires:  python3

%global _description %{expand:
Maximally-Localised Generalised Wannier Functions Code.}

%description
%{_description}

%package        devel
Summary:        Development files for wannier90
Requires:       wannier90%{?_isa} = %{version}-%{release}

%description    devel
This package contains the development files for the wannier90 library.

%package openmpi
Summary:        Maximally-Localised Generalised Wannier Functions Code - OpenMPI version
BuildRequires:  openmpi-devel

%description openmpi
%{_description}

This package contains the OpenMPI parallel version.

%package openmpi-devel
Summary:        Development files for wannier90 - OpenMPI version
Requires:       wannier90-openmpi%{?_isa} = %{version}-%{release}

%description openmpi-devel
This package contains the development files for the wannier90 (OpenMPI) library.

%package mpich
Summary:        Maximally-Localised Generalised Wannier Functions Code - MPICH version
BuildRequires:  mpich-devel

%description mpich
%{_description}

This package contains the MPICH parallel version.

%package mpich-devel
Summary:        Development files for wannier90 - MPICH version
Requires:       wannier90-mpich%{?_isa} = %{version}-%{release}

%description mpich-devel
This package contains the development files for the wannier90 (MPICH) library.


%prep
%autosetup -n wannier90-%{version}

# $MPI_SUFFIX will be evaluated in the loops below, set by mpi modules
%global _vpath_builddir %{_vendor}-%{_target_os}-build${MPI_SUFFIX:-_serial}
# We are running the module load/unload manually until there is a macro-like way to expand this
. /etc/profile.d/modules.sh


%build
cmake_common_args=(
  "-G Ninja"
  "-DWANNIER90_SHARED_LIBS=ON"
  "-DWANNIER90_TEST=ON"
)
for mpi in '' mpich openmpi ; do
  if [ -n "$mpi" ]; then
    module load mpi/${mpi}-%{_arch}
    cmake_mpi_args=(
      "-DCMAKE_INSTALL_PREFIX=${MPI_HOME}"
      "-DWANNIER90_MPI=ON"
      "-DCMAKE_INSTALL_MODULEDIR=${MPI_FORTRAN_MOD_DIR}"
      "-DCMAKE_INSTALL_LIBDIR=lib"
    )
  else
    cmake_mpi_args=(
      "-DWANNIER90_MPI=OFF"
      "-DCMAKE_INSTALL_MODULEDIR=%{_fmoddir}"
    )
  fi

  %cmake \
    ${cmake_common_args[@]} \
    ${cmake_mpi_args[@]}
  %cmake_build

  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done


%install
for mpi in '' mpich openmpi ; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
  %cmake_install
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done


%check
for mpi in '' mpich %{?with_openmpi:openmpi} ; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
  # TODO: re-enable tests
  # %%ctest
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done


%files
%doc README.rst
%license LICENSE
%{_libdir}/libwannier90.so.*
%{_bindir}/wannier90.x
%{_bindir}/postw90.x

%files devel
%{_includedir}/wannier90.hh
%{_libdir}/libwannier90.so
%{_fmoddir}/Wannier90/
%{_libdir}/cmake/Wannier90
%{_libdir}/pkgconfig/wannier90.pc

%files openmpi
%{_libdir}/openmpi/bin/wannier90.x
%{_libdir}/openmpi/bin/postw90.x
%{_libdir}/openmpi/lib/libwannier90.so.*

%files openmpi-devel
%{_libdir}/openmpi/include/wannier90.hh
%{_libdir}/openmpi/lib/libwannier90.so
%{_fmoddir}/openmpi/Wannier90/
%{_libdir}/openmpi/lib/cmake/Wannier90
%{_libdir}/openmpi/lib/pkgconfig/wannier90.pc

%files mpich
%{_libdir}/mpich/bin/wannier90.x
%{_libdir}/mpich/bin/postw90.x
%{_libdir}/mpich/lib/libwannier90.so.*

%files mpich-devel
%{_libdir}/mpich/include/wannier90.hh
%{_libdir}/mpich/lib/libwannier90.so
%{_fmoddir}/mpich/Wannier90/
%{_libdir}/mpich/lib/cmake/Wannier90
%{_libdir}/mpich/lib/pkgconfig/wannier90.pc

%changelog
%autochangelog
