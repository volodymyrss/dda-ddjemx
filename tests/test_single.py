import dataanalysis.core as da

def test_image():
    #
    da.reset()
    import ddosa
    import ddjemx
    reload(ddosa)
    reload(ddjemx)

    da.debug_output()

    fa = ddjemx.jemx_image(assume=[
        ddosa.ScWData(input_scwid="010200230010.001"),
    ])
    fa.read_caches = []

    fa.get()

    print(fa.skyima)

def test_image_nodata():
    #
    da.reset()
    import ddosa
    import ddjemx
    reload(ddosa)
    reload(ddjemx)

    da.debug_output()

    fa = ddjemx.jemx_image(assume=[
        ddosa.ScWData(input_scwid="066500230010.001"),
    ])
    fa.read_caches = []

    try:
        fa.get()
    except ddjemx.ExceptionJ_SCW_NO_MINIMUM_DATA as ex:
        assert ex.args['jemx']=="jmx2"

