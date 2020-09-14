**********************
Environment Variables
**********************

By default :meth:`example_filepath` puts files under the user's cache directory. Override this by setting
the `BGEN_READER_CACHE_HOME` environment variable.

By default, :class:`open_bgen` uses all available processors. Override this with the `num_threads`
parameter or by setting the `MKL_NUM_THREADS` environment variable.
