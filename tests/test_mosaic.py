import dataanalysis.core as da


def test_image_list():
    #da.reset()
    import ddosa
    import ddjemx
    #reload(da)
    reload(ddosa)
    reload(ddjemx)

    jb=ddjemx.JEnergyBins(use_bins=[(3, 10), (10, 30)])
    jb.promote()

    mosaic=ddjemx.JMXScWImageList(
        input_scwlist=ddosa.IDScWList(use_scwid_list=["010200230010.001","010200240010.001","010200250010.001"]),
    )

    mosaic.get()


def test_mosaic():
  #  da.reset()
    import ddosa
    import ddjemx
    #reload(da)
    reload(ddosa)
    reload(ddjemx)


    jb = ddjemx.JEnergyBins(use_bins=[(3, 10), (10, 30)])
    jb.promote()

    mosaic=ddjemx.mosaic_jemx(
              assume=[
                  ddjemx.JMXScWImageList(input_scwlist=ddosa.IDScWList(use_scwid_list=["010200230010.001","010200240010.001","010200250010.001"])),
              ]
            )

    mosaic.get()


def test_mosaic_srcloc():
    #  da.reset()
    import ddosa
    import ddjemx
    # reload(da)
    reload(ddosa)
    reload(ddjemx)

    jb = ddjemx.JEnergyBins(use_bins=[(3, 10), (10, 30)])
    jb.promote()

    src_locator = ddjemx.mosaic_src_loc(input_mosaic=ddjemx.mosaic_jemx(
        assume=[
            ddjemx.JMXScWImageList(input_scwlist=ddosa.IDScWList(
                use_scwid_list=["010200230010.001", "010200240010.001", "010200250010.001"])),
        ]
    ))

    src_locator.get()


def test_spectra_groups():
    #  da.reset()
    import ddosa
    import ddjemx
    # reload(da)
    reload(ddosa)
    reload(ddjemx)

    jb = ddjemx.JEnergyBins(use_bins=[(3, 10), (10, 30)])
    jb.promote()

    groups = ddjemx.JMXImageSpectraGroups(input_scwlist=ddosa.IDScWList(
                use_scwid_list=["010200230010.001", "010200240010.001", "010200250010.001"])
    )

    groups.get()

    groups.construct_og("ogg.fits")

    dl = ddosa.heatool("dal_list")
    dl['dol'] = "ogg.fits"
    dl.run()


def test_spectra_grouped():
    #  da.reset()
    import ddosa
    import ddjemx
    # reload(da)
    reload(ddosa)
    reload(ddjemx)

    jb = ddjemx.JEnergyBins(use_bins=[(3, 10), (10, 30)])
    jb.promote()


    groups = ddjemx.spe_pick(
                use_source_names=["J053432.0+220052"],
                input_spegroups=ddjemx.JMXImageSpectraGroups(input_scwlist=ddosa.IDScWList(
                    use_scwid_list=["010200230010.001", "010200240010.001", "010200250010.001"])
    ))

    groups.get()



