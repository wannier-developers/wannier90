%global soversion 4

Name:           wannier90
Summary:        Maximally-Localised Generalised Wannier Functions Code
Version:        0.0.0
Release:        %autorelease
License:        LGPL-2.1-or-later
URL:            https://www.wannier.org/

Source:         https://github.com/wannier-developers/wannier90/archive/refs/tags/v%{version}.tar.gz
# ix86 because https://fedoraproject.org/wiki/Changes/EncourageI686LeafRemoval
# s390x because tests fail and unknown if upstream wants to support
#  https://github.com/wannier-developers/wannier90/issues/731
ExcludeArch:    %{ix86} s390x

BuildRequires:  cmake
BuildRequires:  gcc-fortran
BuildRequires:  flexiblas-devel
# Required for testing
BuildRequires:  gcc-c++
BuildRequires:  python3dist(pytest)
BuildRequires:  python3dist(pyyaml)

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


%conf
cmake_common_args=(
  "-DWANNIER90_SHARED_LIBS:BOOL=ON"
  "-DWANNIER90_TEST:BOOL=ON"
  "-DWANNIER90_WITH_C:BOOL=ON"
)
for mpi in '' mpich openmpi ; do
  if [ -n "$mpi" ]; then
    module load mpi/${mpi}-%{_arch}
    cmake_mpi_args=(
      "-DCMAKE_INSTALL_PREFIX:PATH=${MPI_HOME}"
      "-DWANNIER90_MPI:BOOL=ON"
      "-DCMAKE_INSTALL_MODULEDIR:PATH=${MPI_FORTRAN_MOD_DIR}"
      "-DCMAKE_INSTALL_LIBDIR:PATH=lib"
    )
  else
    cmake_mpi_args=(
      "-DWANNIER90_MPI:BOOL=OFF"
      "-DCMAKE_INSTALL_MODULEDIR:PATH=%{_fmoddir}"
    )
  fi

  %cmake \
    ${cmake_common_args[@]} \
    ${cmake_mpi_args[@]}

  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done


%build
for mpi in '' mpich openmpi ; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
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
# Skipping a few know test failures for now, tracked in
# https://github.com/wannier-developers/wannier90/issues/666
# https://github.com/wannier-developers/wannier90/issues/731
for mpi in '' mpich openmpi ; do
  [ -n "$mpi" ] && module load mpi/${mpi}-%{_arch}
  %ctest \
    -E "^(library-mode-test-C-interface|testw90_example11_2|testw90_nnkpt4|testw90_nnkpt5)$"
  [ -n "$mpi" ] && module unload mpi/${mpi}-%{_arch}
done


%files
%doc README.rst
%license LICENSE
%{_libdir}/libwannier90.so.%{soversion}{,.*}
%{_bindir}/wannier90.x
%{_bindir}/postw90.x

%files devel
%{_includedir}/wannier90.h
%{_libdir}/libwannier90.so
%{_fmoddir}/Wannier90/
%{_libdir}/cmake/Wannier90
%{_libdir}/pkgconfig/wannier90.pc

%files openmpi
%{_libdir}/openmpi/bin/wannier90.x
%{_libdir}/openmpi/bin/postw90.x
%{_libdir}/openmpi/lib/libwannier90.so.%{soversion}{,.*}

%files openmpi-devel
%{_libdir}/openmpi/include/wannier90.h
%{_libdir}/openmpi/lib/libwannier90.so
%{_fmoddir}/openmpi/Wannier90/
%{_libdir}/openmpi/lib/cmake/Wannier90
%{_libdir}/openmpi/lib/pkgconfig/wannier90.pc

%files mpich
%{_libdir}/mpich/bin/wannier90.x
%{_libdir}/mpich/bin/postw90.x
%{_libdir}/mpich/lib/libwannier90.so.%{soversion}{,.*}

%files mpich-devel
%{_libdir}/mpich/include/wannier90.h
%{_libdir}/mpich/lib/libwannier90.so
%{_fmoddir}/mpich/Wannier90/
%{_libdir}/mpich/lib/cmake/Wannier90
%{_libdir}/mpich/lib/pkgconfig/wannier90.pc

%changelog
%autochangelog
