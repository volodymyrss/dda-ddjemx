import dataanalysis.core as da


def test_image_list():
    da.reset()
    import ddosa
    import ddjemx
    reload(ddosa)
    reload(ddjemx)

    mosaic=ddjemx.JMXScWImageList(input_scwlist=ddosa.IDScWList(use_scwid_list=["010200230010.001","010200240010.001"]))

    mosaic.get()


def test_mosaic():
    da.reset()
    import ddosa
    import ddjemx
    reload(ddosa)
    reload(ddjemx)

    mosaic=ddjemx.mosaic_jemx(
              assume=[
                  ddjemx.JMXScWImageList(input_scwlist=ddosa.IDScWList(use_scwid_list=["010200230010.001","010200240010.001"])),
              ]
            )

    mosaic.get()

