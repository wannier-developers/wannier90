# Library-mode projections regression test

This test exercises the projections parsing path in library mode through `w90_set_option`.

It verifies that when two projection entries are provided, both are parsed correctly.
A previous off-by-one bug in `w90_readwrite_get_projections` skipped the last entry and caused
`w90_input_setopt` to return an error.
