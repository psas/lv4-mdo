# Run Simulation

Run OpenRocket at least once. Use the provided `psas_rocket.ork` file for the
rocket. Ignore the missing motor error.


## Run engine script

Open the `fake-liquid-engine.py` file with a plain text editor and near the top
edit the variables. Save the file with your numbers.

On windows [read about running python programs on windows](https://docs.python.org/2/faq/windows.html)
and run `fake-liquid-engine.py`

This should make a file that OpenRocket can read, so now open, or close then
re-open OpenRocket and the new motor should appear in the select motor dialog.
