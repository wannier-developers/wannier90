#-    if (ierr /= 0) call io_error('Error allocating in_data in w90_readwrite_in_file', stdout, seedname)
#+    if (ierr /= 0) then
#+      call set_error_alloc(error, 'Error allocating in_data in w90_readwrite_in_file')
#+      return
#+    endif

# saves some horse work processing alloc and dealloc failure checks
# takes filename as argument

sed -i "s/\( *\)if (ierr \/= 0).*call io_error(\('.* allocating.*'\), stdout, seedname)/\1if (ierr \/= 0) then\n\1  call set_error_alloc(error, \2)\n\1  return\n\1endif/; s/\( *\)if (ierr \/= 0).*call io_error(\('.* deallocating.*'\), stdout, seedname)/\1if (ierr \/= 0) then\n\1  call set_error_dealloc(error, \2)\n\1  return\n\1endif/" $1
