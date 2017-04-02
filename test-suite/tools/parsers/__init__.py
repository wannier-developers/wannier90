import os

# show_output: "Global" variable that decides whether the parsers should 
# be verbose. Can be changed by setting to 'true' the environment variable
# W90VERBOSETESTS

# Default value
show_output = False

if os.environ.get('W90VERBOSETESTS', 'false').lower() == 'true':
    show_output = True

