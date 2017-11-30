def test_image():
    #
    da.reset()
    import ddosa
    reload(ddosa)

    da.debug_output()

    fa = ddosa.ii_skyimage(assume=[
        ddosa.ScWData(input_scwid="066500230010.001"),
    ])
    fa.read_caches = []

    fa.get()

    print(fa.skyima)