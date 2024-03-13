# %% [markdown]
# 
# 
# ## Summary
# 
# basically a mutilated corpse of the "`swifttools.ukssdc.data.GRB` module" notebook provided by the University of Leister XRT data website when downloading XRT products in bulk/using the swifttools package. It downloads every burst's burst analyzer and could be used to find their plateau start and end times, as well as their afterglow lightcurves' average slope. But everything except downloading the burstanalyzers should be commented out at the moment.

# %%
import swifttools.ukssdc.data.GRB as udg
import os
directory_path = os.getcwd()

# %% [markdown]
# <a id='intro'></a>
# ## Introduction and common behaviour
# 
# Before we dive into things, let me introduce a couple of important concepts which are common to everything that follows.
# 
# Firstly, you specify which GRB(s) you are interested in either by their name, or their targetID. You do this (prepare for a shock) by supplying the `GRBName` or `targetID` argument to the relevant function: all functions in this module to get data take these arguments. You must supply one or the other though, if you suply a `GRBName` and `targetID` you will get an error. These arguments can be either a single value (e.g. `GRBName="GRB 060729"`) or a list/tuple (e.g. `targetID=[282445, 221755]`) if you want more than one GRB. Hold that thought for just a moment because we'll come back to it in a tick.
# 
# The second common concept is how the data you request are stored. This is controlled by two arguments, which again are common to all of the data-access functions in this module. These are:
# 
# * `returnData` - (default `False`), whether or not the function should return a `dict` containing the data.
# * `saveData` - (default: `True`), whether or not the function should download the data files from the UKSSDC website and write them to disk.
# 
# You can set both to `True` to return and save. (You can also set both to `False` if you like wasting time, CPU and making TCP packets feel that they serve no real purpose in life, but why would you do that?) If you are using `saveData=True` you may also need the `clobber` parameter. This specifies whether existing files should be overwritten, and is `False` by default.
# 
# You should still be holding onto a thought&hellip; the fact that you can supply either a single value or a list (or tuple) to identify the GRBs you want. What you do here affects how the data are stored as well.
# 
# If you gave a single value, then `saveData=True` will cause the data for your object to be downloaded to the specified directory. `returnData=True` will cause the function to return a `dict` containing your data.
# 
# If you supplied a list of objects then of course you're going to get multiple dataset. `saveData=True` will (by default) create a subdirectory per GRB to stick the data in. `returnData=True` will return a `dict` with an entry per GRB; that entry will the `dict` containing the data.
# 
# But, you may ask, what will the names of these subdirectories, or the keys of these `dict`s be? They will be whatever identifer you used to specify the GRBs you wanted. If you used the `GRBName` argument, the directories and `dict` keys will be the GRB names; if you used `targetID` then they will be, er, the target IDs.
# 
# That may not make so much sense in the abstract, but all will become clear as we see some demos. Oh, but just to further confuse you, please do note that a list (or tuple) with one entry is still a list (unless it's a tuple), and so the returned data will still have the subdirectories/extra `dict` that you get when you supply a list (or tuple) even though there is only one entry in this. Confused? Me too, but it will make sense when we start playing, so let's do so.

# %%
lcData = udg.getLightCurves(GRBName=("GRB 230414B", "GRB230409B", "GRB230405B",\
                                              "GRB230328B", "GRB230325A", "GRB230322B",\
                                              "GRB230307A", "GRB230228A", "GRB230217A",\
                                              "GRB230216A", "GRB230205A", \
                                              "GRB230204B", "GRB230204A", "GRB230123A",\
                                              "GRB230116D", "GRB221231A", "GRB221226B",\
                                              "GRB221216A", "GRB221202A", "GRB221201A",\
                                              "GRB221120A", "GRB221110A", "GRB221028A", \
                                              "GRB221027B", "GRB221024A", "GRB221016A", \
                                              "GRB221009A", "GRB221006A", "GRB220930A",\
                                              "GRB220921A", "GRB220907A", "GRB220826A",\
                                              "GRB220813A", "GRB220730A", "GRB220715B",\
                                              "GRB220714B", "GRB220711B", "GRB220710A",\
                                              "GRB220708B", "GRB220708A", "GRB220706A", \
                                              "GRB220701A", "GRB220627A", "GRB220623A", \
                                              "GRB220618A", "GRB220611A", "GRB220527A",\
                                              "GRB220521A", "GRB220518A", "GRB220511A", \
                                              "GRB220506A", "GRB220430A", "GRB220427A", \
                                              "GRB220412A", "GRB220408A", "GRB220403B", \
                                              "GRB220325A", "GRB220319A", "GRB220306B",\
                                              "GRB220305A", "GRB220302A", "GRB220219B",\
                                              "GRB220118A", "GRB220117B", "GRB220117A", \
                                              "GRB220107B", "GRB220107A", "GRB220101A",\
                                              "GRB211227A", "GRB211221A", "GRB211211A", \
                                              "GRB211207A", "GRB211129A", "GRB211107B",\
                                              "GRB211106A", "GRB211025A", "GRB211024B",\
                                              "GRB211023B", "GRB210930A", "GRB210928A",\
                                              "GRB210919A", "GRB210912A", "GRB210905A", \
                                              "GRB210901A", "GRB210827A", "GRB210824A",\
                                              "GRB210822A", "GRB210820A", "GRB210818A",\
                                              "GRB210807C", "GRB210807A", "GRB210802A", \
                                              "GRB210731A", "GRB210730A", "GRB210726A", \
                                              "GRB210725B", "GRB210725A", "GRB210724A", \
                                              "GRB210723A", "GRB210722A", "GRB210712A", \
                                              "GRB210708A", "GRB210706A", "GRB210704A", \
                                              "GRB210702A", "GRB210626A", "GRB210619B", \
                                              "GRB210618A", "GRB210610B", "GRB210610A",\
                                              "GRB210606A", "GRB210605B", "GRB210527A", \
                                              "GRB210517A", "GRB210515C", "GRB210514A", \
                                              "GRB210509A", "GRB210504A", "GRB210421A", \
                                              "GRB210420B", "GRB210419C", "GRB210419A", \
                                              "GRB210411C", "GRB210410A", "GRB210402A", \
                                              "GRB210323A", "GRB210321A", "GRB210318B", \
                                              "GRB210318A", "GRB210308A", "GRB210307A", \
                                              "GRB210306A", "GRB210305A", "GRB210226A", \
                                              "GRB210222B", "GRB210218A", "GRB210217A", \
                                              "GRB210212A", "GRB210211A", "GRB210210A", \
                                              "GRB210209A", "GRB210207B", "GRB210205A", \
                                              "GRB210204A", "GRB210112A", "GRB210104B", \
                                              "GRB210104A", "GRB210102B", "GRB201229A", \
                                              "GRB201223A", "GRB201221D", "GRB201221A", \
                                              "GRB201216C", "GRB201209A", "GRB201208A", \
                                              "GRB201203A", "GRB201128A", "GRB201116A", \
                                              "GRB201104B", "GRB201103B", "GRB201029A",\
                                              "GRB201027A", "GRB201026A", "GRB201024A", \
                                              "GRB201021C", "GRB201020B", "GRB201020A", \
                                              "GRB201017A", "GRB201015A", "GRB201014A", \
                                              "GRB201013A", "GRB201008A", "GRB201006A", \
                                              "GRB201001A", "GRB200925B", "GRB200922A",\
                                              "GRB200917A", "GRB200907B", "GRB200906A", \
                                              "GRB200901B", "GRB200901A", "GRB200829A", \
                                              "GRB200826A", "GRB200819A", "GRB200809B", \
                                              "GRB200806A", "GRB200803A", "GRB200729A", \
                                              "GRB200716C", "GRB200715A", "GRB200714E", \
                                              "GRB200713A", "GRB200711A", "GRB200630A", \
                                              "GRB200613A", "GRB200612A", "GRB200529A",\
                                              "GRB200528A", "GRB200524A", "GRB200522A", \
                                              "GRB200519A", "GRB200517A", "GRB200512A", \
                                              "GRB200509A", "GRB200425A", "GRB200424A", \
                                              "GRB200416A", "GRB200412B", "GRB200411A",\
                                              "GRB200410A", "GRB200409A", "GRB200324A", \
                                              "GRB200306C", "GRB200306A", "GRB200303A", \
                                              "GRB200228B", "GRB200227A", "GRB200224A", \
                                              "GRB200219C", "GRB200219A", "GRB200216B", \
                                              "GRB200215A", "GRB200205B", "GRB200205A", \
                                              "GRB200131A", "GRB200127A", "GRB200125A", \
                                              "GRB200122A", "GRB200109A", "GRB200107B", \
                                              "GRB191228A", "GRB191221B", "GRB191220A", \
                                              "GRB191123A", "GRB191122A", "GRB191106A", \
                                              "GRB191101A", "GRB191031D", "GRB191031C", \
                                              "GRB191029A", "GRB191024A", "GRB191019A", \
                                              "GRB191017B", "GRB191016A", "GRB191011A", \
                                              "GRB191004B", "GRB191004A", "GRB190926A", \
                                              "GRB190919B", "GRB190829A", "GRB190828B", \
                                              "GRB190824A", "GRB190823A", "GRB190821A", \
                                              "GRB190816A", "GRB190804B", "GRB190731A", \
                                              "GRB190719C", "GRB190718A", "GRB190708B", \
                                              "GRB190706B", "GRB190701A", "GRB190630C", \
                                              "GRB190630B", "GRB190627A", "GRB190613B", \
                                              "GRB190613A", "GRB190611A", "GRB190604B", \
                                              "GRB190531B", "GRB190530A", "GRB190519A", \
                                              "GRB190515B", "GRB190511A", "GRB190422A", \
                                              "GRB190404C", "GRB190324A", "GRB190320A", \
                                              "GRB190311A", "GRB190305A", "GRB190220A", \
                                              "GRB190219A", "GRB190211A", "GRB190204A", \
                                              "GRB190203A", "GRB190202A", "GRB190129B", \
                                              "GRB190123A", "GRB190114C", "GRB190114B", \
                                              "GRB190114A", "GRB190109B", "GRB190109A", \
                                              "GRB190106A", "GRB190103B", "GRB181228A", \
                                              "GRB181213A", "GRB181203A", "GRB181202A", \
                                              "GRB181201A", "GRB181126A", "GRB181123B", \
                                              "GRB181110A", "GRB181103A", "GRB181030A", \
                                              "GRB181023A", "GRB181022A", "GRB181020A", \
                                              "GRB181013A", "GRB181010A", "GRB181003A", \
                                              "GRB181002A", "GRB180930A", "GRB180925A", \
                                              "GRB180924A", "GRB180914B", "GRB180905A", \
                                              "GRB180904A", "GRB180828A", "GRB180823A", \
                                              "GRB180821A", "GRB180818B", "GRB180818A", \
                                              "GRB180812A", "GRB180809B", "GRB180809A", \
                                              "GRB180806A", "GRB180805B", "GRB180805A",\
                                              "GRB180728A", "GRB180727A", "GRB180721A", \
                                              "GRB180720C", "GRB180720B", "GRB180709A", \
                                              "GRB180706A", "GRB180705A", "GRB180704A", \
                                              "GRB180703A", "GRB180630A", "GRB180626A", \
                                              "GRB180624A", "GRB180623A", "GRB180620B", \
                                              "GRB180620A", "GRB180618A", "GRB180614A", \
                                              "GRB180613A", "GRB180602A", "GRB180514A", \
                                              "GRB180512A", "GRB180510B", "GRB180510A", \
                                              "GRB180504A", "GRB180425A", "GRB180418A", \
                                              "GRB180411A", "GRB180410A", "GRB180409A", \
                                              "GRB180404B", "GRB180404A", "GRB180402A", \
                                              "GRB180331B", "GRB180331A", "GRB180329B", \
                                              "GRB180325A", "GRB180324A", "GRB180316A", \
                                              "GRB180314B", "GRB180314A", "GRB180311A", \
                                              "GRB180305A", "GRB180224A", "GRB180222A", \
                                              "GRB180210A", "GRB180205A", "GRB180204A", \
                                              "GRB180115A", "GRB180111A", "GRB180103A", \
                                              "GRB180102A", "GRB171222A", "GRB171216A", \
                                              "GRB171212A", "GRB171211A", "GRB171209A", \
                                              "GRB171205A", "GRB171124A", "GRB171123A", \
                                              "GRB171120A", "GRB171115A", "GRB171102B", \
                                              "GRB171027A", "GRB171020A", "GRB171010B", \
                                              "GRB171010A", "GRB171007A", "GRB171004A", \
                                              "GRB171001A", "GRB170921A", "GRB170912B", \
                                              "GRB170912A", "GRB170906C", "GRB170906B", \
                                              "GRB170906A", "GRB170903A", "GRB170827A", \
                                              "GRB170822A", "GRB170813A", "GRB170810A", \
                                              "GRB170804A", "GRB170803A", "GRB170728B", \
                                              "GRB170728A", "GRB170714A", "GRB170711A", \
                                              "GRB170710A", "GRB170705A", "GRB170629A", \
                                              "GRB170626A", "GRB170607A", "GRB170604A", \
                                              "GRB170531B", "GRB170531A", "GRB170526A", \
                                              "GRB170524B", "GRB170524A", "GRB170519A", \
                                              "GRB170516A", "GRB170510A", "GRB170428A", \
                                              "GRB170419A", "GRB170405A", "GRB170331A", \
                                              "GRB170330A", "GRB170318B", "GRB170318A", \
                                              "GRB170317A", "GRB170311A", "GRB170306A", \
                                              "GRB170214A", "GRB170208B", "GRB170208A", \
                                              "GRB170206B", "GRB170205A", "GRB170202A", \
                                              "GRB170127B", "GRB170127A", "GRB170126A", \
                                              "GRB170115A", "GRB170113A", "GRB170111A", \
                                              "GRB161224A", "GRB161220A", "GRB161219B", \
                                              "GRB161217A", "GRB161214B", "GRB161202A", \
                                              "GRB161129A", "GRB161117B", "GRB161117A", \
                                              "GRB161113A", "GRB161108A", "GRB161105A", \
                                              "GRB161104A", "GRB161023A", "GRB161022A", \
                                              "GRB161017A", "GRB161015A", "GRB161014A", \
                                              "GRB161011A", "GRB161010A", "GRB161007A", \
                                              "GRB161004B", "GRB161004A", "GRB161001A", \
                                              "GRB160927A", "GRB160917A", "GRB160912A", \
                                              "GRB160910A", "GRB160905A", "GRB160826A", \
                                              "GRB160824A", "GRB160821B", "GRB160816A", \
                                              "GRB160815A", "GRB160804A", "GRB160801A", \
                                              "GRB160716A", "GRB160712A", "GRB160705B", \
                                              "GRB160703A", "GRB160630A", "GRB160629A", \
                                              "GRB160625B", "GRB160625A", "GRB160624A", \
                                              "GRB160623A", "GRB160611A", "GRB160607A", \
                                              "GRB160601A", "GRB160525B", "GRB160521B", \
                                              "GRB160509A", "GRB160506A", "GRB160504A", \
                                              "GRB160501A", "GRB160425A", "GRB160422A", \
                                              "GRB160417A", "GRB160412A", "GRB160411A", \
                                              "GRB160410A", "GRB160408A", "GRB160327A", \
                                              "GRB160325A", "GRB160321A", "GRB160314A", \
                                              "GRB160313A", "GRB160303A", "GRB160228A", \
                                              "GRB160227A", "GRB160225A", "GRB160223A", \
                                              "GRB160221A", "GRB160220B", "GRB160220A", \
                                              "GRB160216A", "GRB160203A", "GRB160131A", \
                                              "GRB160127A", "GRB160123A", "GRB160121A", \
                                              "GRB160119A", "GRB160117B", "GRB160117A", \
                                              "GRB160104A", "GRB160101A", "GRB151229A", \
                                              "GRB151228B", "GRB151215A", "GRB151212A", \
                                              "GRB151210A", "GRB151205A", "GRB151127A", \
                                              "GRB151120A", "GRB151118A", "GRB151114A", \
                                              "GRB151112A", "GRB151111A", "GRB151031A", \
                                              "GRB151029A", "GRB151027B", "GRB151027A", \
                                              "GRB151023A", "GRB151022A", "GRB151021A", \
                                              "GRB151006A", "GRB151004A", "GRB151001B", \
                                              "GRB151001A", "GRB150925A", "GRB150915A", \
                                              "GRB150911A", "GRB150910A", "GRB150907B", \
                                              "GRB150902A", "GRB150831B", "GRB150831A", \
                                              "GRB150821A", "GRB150819A", "GRB150818A", \
                                              "GRB150817A", "GRB150811A", "GRB150801B", \
                                              "GRB150728A", "GRB150727A", "GRB150724B", \
                                              "GRB150724A", "GRB150722A", "GRB150720A", \
                                              "GRB150716A", "GRB150711A", "GRB150710B", \
                                              "GRB150627A", "GRB150626B", "GRB150626A", \
                                              "GRB150622A", "GRB150616A", "GRB150615A", \
                                              "GRB150608A", "GRB150607A", "GRB150530A", \
                                              "GRB150527A", "GRB150523A", "GRB150518A", \
                                              "GRB150514A", "GRB150430A", "GRB150428B", \
                                              "GRB150428A", "GRB150424A", "GRB150423A", \
                                              "GRB150403A", "GRB150323C", "GRB150323B", \
                                              "GRB150323A", "GRB150318A", "GRB150317A", \
                                              "GRB150314A", "GRB150309A", "GRB150302A", \
                                              "GRB150301B", "GRB150301A", "GRB150222A", \
                                              "GRB150219A", "GRB150213B", "GRB150212A", \
                                              "GRB150211A", "GRB150206A", "GRB150204A", \
                                              "GRB150203A", "GRB150202B", "GRB150202A", \
                                              "GRB150201A", "GRB150120B", "GRB150120A", \
                                              "GRB150110B", "GRB150103A", "GRB150101B", \
                                              "GRB150101A", "GRB141225A", "GRB141221A", \
                                              "GRB141220A", "GRB141215A", "GRB141212B", \
                                              "GRB141212A", "GRB141130A", "GRB141121A", \
                                              "GRB141109B", "GRB141109A", "GRB141031B", \
                                              "GRB141031A", "GRB141028A", "GRB141026A", \
                                              "GRB141022A", "GRB141020A", "GRB141017A", \
                                              "GRB141015A", "GRB141005A", "GRB141004A", \
                                              "GRB140930B", "GRB140928A", "GRB140927A", \
                                              "GRB140919A", "GRB140916A", "GRB140909A", \
                                              "GRB140907A", "GRB140903A", "GRB140824A", \
                                              "GRB140818B", "GRB140818A", "GRB140817A", \
                                              "GRB140815A", "GRB140808A", "GRB140801A", \
                                              "GRB140730A", "GRB140719A", "GRB140716A", \
                                              "GRB140713A", "GRB140710A", "GRB140709B", \
                                              "GRB140709A", "GRB140706A", "GRB140703A", \
                                              "GRB140629A", "GRB140628A", "GRB140626A", \
                                              "GRB140623A", "GRB140622A", "GRB140620A", \
                                              "GRB140619B", "GRB140619A", "GRB140614B", \
                                              "GRB140614A", "GRB140611A", "GRB140610A", \
                                              "GRB140606B", "GRB140529A", "GRB140518A", \
                                              "GRB140516A", "GRB140515A", "GRB140512A", \
                                              "GRB140509A", "GRB140508A", "GRB140506A", \
                                              "GRB140502A", "GRB140430A", "GRB140428A", \
                                              "GRB140423A", "GRB140419A", "GRB140413A", \
                                              "GRB140412A", "GRB140408A", "GRB140331A", \
                                              "GRB140323A", "GRB140320C", "GRB140320B", \
                                              "GRB140320A", "GRB140318A", "GRB140311B", \
                                              "GRB140311A", "GRB140304A", "GRB140302A", \
                                              "GRB140301A", "GRB140226A", "GRB140215A", \
                                              "GRB140213A", "GRB140211A", "GRB140209A", \
                                              "GRB140206A", "GRB140129B", "GRB140129A", \
                                              "GRB140114A", "GRB140108A", "GRB140104B", \
                                              "GRB140103A", "GRB140102A", "GRB131231A", \
                                              "GRB131229A", "GRB131227A", "GRB131205A", \
                                              "GRB131202A", "GRB131128A", "GRB131127A", \
                                              "GRB131122A", "GRB131117A", "GRB131108A", \
                                              "GRB131105A", "GRB131103A", "GRB131031A", \
                                              "GRB131030A", "GRB131024B", "GRB131024A", \
                                              "GRB131018B", "GRB131018A", "GRB131014A", \
                                              "GRB131011A", "GRB131004A", "GRB131002B", \
                                              "GRB131002A", "GRB130929A", "GRB130925A", \
                                              "GRB130912A", "GRB130907A", "GRB130903A", \
                                              "GRB130831B", "GRB130831A", "GRB130822A", \
                                              "GRB130816B", "GRB130816A", "GRB130812A", \
                                              "GRB130807A", "GRB130806A", "GRB130803A", \
                                              "GRB130727A", "GRB130725B", "GRB130725A", \
                                              "GRB130722A", "GRB130716A", "GRB130702A", \
                                              "GRB130701A", "GRB130627B", "GRB130627A", \
                                              "GRB130625A", "GRB130623A", "GRB130615A", \
                                              "GRB130612A", "GRB130610A", "GRB130609B", \
                                              "GRB130609A", "GRB130608A", "GRB130606B", \
                                              "GRB130606A", "GRB130605A", "GRB130604A", \
                                              "GRB130603B", "GRB130603A", "GRB130529A", \
                                              "GRB130528A", "GRB130527A", "GRB130518A", \
                                              "GRB130515A", "GRB130514B", "GRB130514A", \
                                              "GRB130513A", "GRB130511A", "GRB130508A", \
                                              "GRB130505B", "GRB130505A", "GRB130504C", \
                                              "GRB130504A", "GRB130502B", "GRB130502A", \
                                              "GRB130427B", "GRB130427A", "GRB130420B", \
                                              "GRB130420A", "GRB130418A", "GRB130408A", \
                                              "GRB130327B", "GRB130327A", "GRB130315A", \
                                              "GRB130313A", "GRB130306A", "GRB130305A", \
                                              "GRB130211A", "GRB130206A", "GRB130131B", \
                                              "GRB130131A", "GRB130122A", "GRB130102A", \
                                              "GRB121229A", "GRB121226A", "GRB121217A", \
                                              "GRB121212A", "GRB121211A", "GRB121209A", \
                                              "GRB121202A", "GRB121201A", "GRB121128A", \
                                              "GRB121125A", "GRB121123A", "GRB121117A", \
                                              "GRB121108A", "GRB121102A", "GRB121031A", \
                                              "GRB121028A", "GRB121027A", "GRB121025A", \
                                              "GRB121024A", "GRB121017A", "GRB121011A", \
                                              "GRB121001A", "GRB120927A", "GRB120923A", \
                                              "GRB120922A", "GRB120911A", "GRB120909A", \
                                              "GRB120907A", "GRB120819A", "GRB120817A", \
                                              "GRB120816A", "GRB120815A", "GRB120811C", \
                                              "GRB120811A", "GRB120807A", "GRB120805A", \
                                              "GRB120804A", "GRB120803B", "GRB120803A", \
                                              "GRB120802A", "GRB120729A", "GRB120728A", \
                                              "GRB120724A", "GRB120722A", "GRB120716A", \
                                              "GRB120714A", "GRB120712A", "GRB120711B", \
                                              "GRB120711A", "GRB120703A", "GRB120701A", \
                                              "GRB120630A", "GRB120624B", "GRB120612A", \
                                              "GRB120521C", "GRB120521B", "GRB120521A", \
                                              "GRB120514A", "GRB120422A", "GRB120419A", \
                                              "GRB120404A", "GRB120403B", "GRB120401A", \
                                              "GRB120328A", "GRB120327A", "GRB120326A", \
                                              "GRB120324A", "GRB120320A", "GRB120312A", \
                                              "GRB120311B", "GRB120311A", "GRB120308A", \
                                              "GRB120305A", "GRB120302A", "GRB120224A", \
                                              "GRB120219A", "GRB120215A", "GRB120213A", \
                                              "GRB120212A", "GRB120211A", "GRB120202A", \
                                              "GRB120121A", "GRB120119A", "GRB120118B", \
                                              "GRB120116A", "GRB120106A", "GRB120102A", \
                                              "GRB111229A", "GRB111228A", "GRB111225A", \
                                              "GRB111222A", "GRB111215B", "GRB111215A", \
                                              "GRB111212A", "GRB111211A", "GRB111210A", \
                                              "GRB111209A", "GRB111208A", "GRB111205A", \
                                              "GRB111204A", "GRB111201A", "GRB111129A", \
                                              "GRB111123A", "GRB111121A", "GRB111117A", \
                                              "GRB111109A", "GRB111107A", "GRB111103B", \
                                              "GRB111029A", "GRB111022B", "GRB111022A", \
                                              "GRB111020A", "GRB111018A", "GRB111016A", \
                                              "GRB111008A", "GRB110928A", "GRB110921A", \
                                              "GRB110918A", "GRB110915B", "GRB110915A", \
                                              "GRB110903A", "GRB110820A", "GRB110818A", \
                                              "GRB110808A", "GRB110802A", "GRB110801A", \
                                              "GRB110731A", "GRB110726A", "GRB110719A", \
                                              "GRB110715A", "GRB110709B", "GRB110709A", \
                                              "GRB110708A", "GRB110625A", "GRB110610A", \
                                              "GRB110604A", "GRB110530A", "GRB110521A", \
                                              "GRB110520A", "GRB110503A", "GRB110428A", \
                                              "GRB110422A", "GRB110420A", "GRB110414A", \
                                              "GRB110411A", "GRB110407A", "GRB110402A", \
                                              "GRB110319B", "GRB110319A", "GRB110318B", \
                                              "GRB110315A", "GRB110312A", "GRB110305A", \
                                              "GRB110223B", "GRB110223A", "GRB110213B", \
                                              "GRB110213A", "GRB110210A", "GRB110208A", \
                                              "GRB110206A", "GRB110205A", "GRB110201A", \
                                              "GRB110128A", "GRB110119A", "GRB110112A", \
                                              "GRB110107A", "GRB110106B", "GRB110106A", \
                                              "GRB110102A", "GRB101225A", "GRB101224A", \
                                              "GRB101219B", "GRB101219A", "GRB101213A", \
                                              "GRB101204A", "GRB101201A", "GRB101117B", \
                                              "GRB101114A", "GRB101112A", "GRB101030A", \
                                              "GRB101024A", "GRB101023A", "GRB101017A", \
                                              "GRB101011A", "GRB101008A", "GRB100915A", \
                                              "GRB100909A", "GRB100906A", "GRB100905A", \
                                              "GRB100902A", "GRB100901A", "GRB100823A", \
                                              "GRB100816A", "GRB100814A", "GRB100807A", \
                                              "GRB100805A", "GRB100802A", "GRB100728B", \
                                              "GRB100728A", "GRB100727A", "GRB100725B", \
                                              "GRB100725A", "GRB100724A", "GRB100713A", \
                                              "GRB100704A", "GRB100702A", "GRB100628A", \
                                              "GRB100625A", "GRB100621A", "GRB100619A", \
                                              "GRB100615A", "GRB100614A", "GRB100606A", \
                                              "GRB100528A", "GRB100526B", "GRB100526A", \
                                              "GRB100522A", "GRB100518A", "GRB100514A", \
                                              "GRB100513A", "GRB100508A", "GRB100504A", \
                                              "GRB100425A", "GRB100424A", "GRB100420A", \
                                              "GRB100418A", "GRB100414A", "GRB100413A", \
                                              "GRB100331B", "GRB100316D", "GRB100316C", \
                                              "GRB100316B", "GRB100316A", "GRB100305A", \
                                              "GRB100302A", "GRB100219A", "GRB100213B", \
                                              "GRB100213A", "GRB100212A", "GRB100206A", \
                                              "GRB100205A", "GRB100117A", "GRB100115A", \
                                              "GRB100111A", "GRB100103A", "GRB091230", \
                                              "GRB091221", "GRB091208B", "GRB091202", \
                                              "GRB091130B", "GRB091127", "GRB091111", \
                                              "GRB091109B", "GRB091109A", "GRB091104", \
                                              "GRB091102", "GRB091029", "GRB091026", \
                                              "GRB091024", "GRB091020", "GRB091018", \
                                              "GRB091010", "GRB091003", "GRB090929B", \
                                              "GRB090927", "GRB090926B", "GRB090926A", \
                                              "GRB090915", "GRB090912", "GRB090904B", \
                                              "GRB090904A", "GRB090902B", "GRB090831C", \
                                              "GRB090827", "GRB090823", "GRB090817", \
                                              "GRB090814B", "GRB090814A", "GRB090813", \
                                              "GRB090812", "GRB090809", "GRB090807", \
                                              "GRB090728", "GRB090727", "GRB090726",\
                                              "GRB090720", "GRB090715B", "GRB090709A", \
                                              "GRB090702", "GRB090628", "GRB090625B", \
                                              "GRB090621B", "GRB090621A", "GRB090618", \
                                              "GRB090607", "GRB090531B", "GRB090531A", \
                                              "GRB090530", "GRB090529", "GRB090519", \
                                              "GRB090518", "GRB090516", "GRB090515", \
                                              "GRB090510", "GRB090429B", "GRB090429A", \
                                              "GRB090426", "GRB090424", "GRB090423", \
                                              "GRB090422", "GRB090419", "GRB090418A", \
                                              "GRB090417B", "GRB090407", "GRB090404", \
                                              "GRB090401B", "GRB090328A", "GRB090323", \
                                              "GRB090313", "GRB090309", "GRB090308", \
                                              "GRB090307", "GRB090306B", "GRB090205", \
                                              "GRB090201", "GRB090126A", "GRB090123", \
                                              "GRB090118", "GRB090117", "GRB090113", \
                                              "GRB090111", "GRB090107B", "GRB090102", \
                                              "GRB081230", "GRB081228", "GRB081226A", \
                                              "GRB081224", "GRB081222", "GRB081221", \
                                              #"GRB081211B", 
                                                 "GRB081211", "GRB081210", \
                                              "GRB081204", "GRB081203B", "GRB081203A", \
                                              "GRB081128", "GRB081127", "GRB081126", \
                                              "GRB081121", "GRB081118", "GRB081109", \
                                              "GRB081105", "GRB081104", "GRB081102", \
                                              "GRB081029", "GRB081028", "GRB081025", \
                                              "GRB081024A", "GRB081016B", "GRB081016A", \
                                              "GRB081012", "GRB081011", "GRB081008", \
                                              "GRB081007", "GRB081003A", "GRB081001", \
                                              "GRB080928", "GRB080919", "GRB080916C", \
                                              "GRB080916B", "GRB080916A", "GRB080915A", \
                                              "GRB080913", "GRB080906", "GRB080905B", \
                                              "GRB080905A", "GRB080903", "GRB080828", \
                                              "GRB080825B", "GRB080810", "GRB080805", \
                                              "GRB080804", "GRB080802", "GRB080727C", \
                                              "GRB080727B", "GRB080727A", "GRB080723B", \
                                              "GRB080723A", "GRB080721", "GRB080714", \
                                              "GRB080710", "GRB080707", "GRB080703", \
                                              "GRB080702B", "GRB080702A", "GRB080701", \
                                              "GRB080625", "GRB080623", "GRB080613B", \
                                              "GRB080613A", "GRB080607", "GRB080605", \
                                              "GRB080604", "GRB080603B", "GRB080603A", \
                                              "GRB080602", "GRB080523", "GRB080520", \
                                              "GRB080517", "GRB080516", "GRB080515", \
                                              "GRB080514B", "GRB080507", "GRB080506", \
                                              "GRB080503", "GRB080430", "GRB080426", \
                                              "GRB080413B", "GRB080413A", "GRB080411", \
                                              "GRB080409", "GRB080405", "GRB080330", \
                                              "GRB080328", "GRB080325", "GRB080320", \
                                              "GRB080319D", "GRB080319C", "GRB080319B", \
                                              "GRB080319A", "GRB080310", "GRB080307", \
                                              "GRB080303", "GRB080229B", "GRB080229A", \
                                              "GRB080218B", "GRB080212", "GRB080210", \
                                              "GRB080207", "GRB080205", "GRB080130", \
                                              "GRB080129", "GRB080123", "GRB080120", \
                                              "GRB071227", "GRB071122", "GRB071118", \
                                              "GRB071117", "GRB071112C", "GRB071104", \
                                              "GRB071101", "GRB071031", "GRB071028B", \
                                              "GRB071028A", "GRB071025", "GRB071021", \
                                              "GRB071020", "GRB071011", "GRB071010B", \
                                              "GRB071010A", "GRB071008", "GRB071003", \
                                              "GRB070925", "GRB070920B", "GRB070917", \
                                              "GRB070913", "GRB070911", "GRB070810A", \
                                              "GRB070809", "GRB070808", "GRB070802", \
                                              "GRB070729", "GRB070724B", "GRB070724A", \
                                              "GRB070721B", "GRB070721A", "GRB070714B", \
                                              "GRB070714A", "GRB070707", "GRB070704", \
                                              "GRB070628", "GRB070621", "GRB070616", \
                                              "GRB070615", "GRB070612B", "GRB070611", \
                                              "GRB070531", "GRB070529", "GRB070521", \
                                              "GRB070520B", "GRB070520A", "GRB070518", \
                                              "GRB070517", "GRB070509", "GRB070508", \
                                              "GRB070506", "GRB070429B", "GRB070429A", \
                                              "GRB070420", "GRB070419B", "GRB070419A", \
                                              "GRB070412", "GRB070411", "GRB070330", \
                                              "GRB070328", "GRB070318", "GRB070311", \
                                              "GRB070309", "GRB070306", "GRB070227", \
                                              "GRB070224", "GRB070223", "GRB070220", \
                                              "GRB070219", "GRB070208", "GRB070129", \
                                              "GRB070125", "GRB070110", "GRB070107", \
                                              "GRB070103", "GRB061222B", "GRB061222A", \
                                              "GRB061217", "GRB061210", "GRB061202", \
                                              "GRB061201", "GRB061126", "GRB061122", \
                                              "GRB061121", "GRB061110B", "GRB061110A", \
                                              "GRB061102", "GRB061028", "GRB061025", \
                                              "GRB061021", "GRB061019", "GRB061007", \
                                              "GRB061006", "GRB061004", "GRB061002", \
                                              "GRB060929", "GRB060928", "GRB060927", \
                                              "GRB060926", "GRB060923C", "GRB060923B", \
                                              "GRB060923A", "GRB060919", "GRB060912A", \
                                              "GRB060908", "GRB060906", "GRB060904B", \
                                              "GRB060904A", "GRB060901", "GRB060825", \
                                              "GRB060814", "GRB060813", "GRB060807", \
                                              "GRB060805B", "GRB060805A", "GRB060804", \
                                              "GRB060801", "GRB060729", "GRB060719", \
                                              "GRB060717", "GRB060714", "GRB060712", \
                                              "GRB060708", "GRB060707", "GRB060614", \
                                              "GRB060607A", "GRB060605", "GRB060604", \
                                              "GRB060602B", "GRB060602A", "GRB060526", \
                                              "GRB060522", "GRB060515", "GRB060512", \
                                              "GRB060510B", "GRB060510A", "GRB060507", \
                                              "GRB060505", "GRB060502B", "GRB060502A", \
                                              "GRB060501", "GRB060428B", "GRB060428A", \
                                              "GRB060427", "GRB060421", "GRB060418", \
                                              "GRB060413", "GRB060403", "GRB060323", \
                                              "GRB060322", "GRB060319", "GRB060313", \
                                              "GRB060312", "GRB060306", "GRB060223A", \
                                              "GRB060219", "GRB060218", "GRB060211B", \
                                              "GRB060211A", "GRB060210", "GRB060206", \
                                              "GRB060204B", "GRB060203", "GRB060202", \
                                              "GRB060124", "GRB060123", "GRB060121", \
                                              "GRB060116", "GRB060115", "GRB060111B", \
                                              "GRB060111A", "GRB060110", "GRB060109", \
                                              "GRB060108", "GRB060105", "GRB051227", \
                                              "GRB051221B", "GRB051221A", "GRB051211B",\
                                              "GRB051210", "GRB051117B", "GRB051117A", \
                                              "GRB051111", "GRB051109B", "GRB051109A", \
                                              "GRB051028", "GRB051022", "GRB051021B", \
                                              "GRB051021A", "GRB051016B", "GRB051016A", \
                                              "GRB051012", "GRB051008", "GRB051006", \
                                              "GRB051001", "GRB050922C", "GRB050922B", \
                                              "GRB050916", "GRB050915B", "GRB050915A", \
                                              "GRB050908", "GRB050904", "GRB050827", \
                                              "GRB050826", "GRB050824", "GRB050822", \
                                              "GRB050820B", "GRB050820A", "GRB050819", \
                                              "GRB050815", "GRB050814", "GRB050813", \
                                              "GRB050803", "GRB050802", "GRB050801", \
                                              "GRB050730", "GRB050726", "GRB050724", \
                                              "GRB050721", "GRB050717", "GRB050716", \
                                              "GRB050714B", "GRB050714", "GRB050713B", \
                                              "GRB050713A", "GRB050712", "GRB050701", \
                                              "GRB050607", "GRB050603", "GRB050525A", \
                                              "GRB050520", "GRB050509B", "GRB050509A", \
                                              "GRB050505", "GRB050504", "GRB050502B", \
                                              "GRB050422", "GRB050421", "GRB050416A", \
                                              "GRB050412", "GRB050410", "GRB050408", \
                                              "GRB050406", "GRB050401", "GRB050326", \
                                              "GRB050319", "GRB050318", "GRB050315", \
                                              "GRB050223", "GRB050219B", "GRB050219A", \
                                              "GRB050215B", "GRB050128", "GRB050126", \
                                              "GRB050124", "GRB041223", "GRB041218"),
                            destDir='{}/xrt_lcs/'('{}/xrt_lcs/{}'.format(directory_path,,
                            subDirs=False,
                            silent=False,
                            verbose=True,
                            saveData=True
                            )

# %% [markdown]
# ### Storing the light curves in variables
# 
# Let's move on now to the `returnData=True` case. As I told you [earlier](#intro) this will return a `dict` containing the data. All light curves returned by anything in the `swifttools.ukssdc` module have a common structure which I call a "light curve `dict`", and you can [read about this here](https://www.swift.ac.uk/API/ukssdc/structures.md#the-light-curve-dict).
# 
# There are no special parameters related to returning data, so let's jump straight in with some demonstrations similar to those above.

# %%
# lcData = udg.getLightCurves(GRBName="GRB 220427A",
#                             incbad="both",
#                             nosys="both",
#                             saveData=False,
#                             returnData=True)

# %% [markdown]
# I don't recommend simply printing `lcData` straight, it's quite big. If you [read about the light curve `dict`](https://www.swift.ac.uk/API/ukssdc/structures.md#the-light-curve-dict) then you may have an idea what to expect, but it's helpful to explore it anyway, so let's do that.

# %%
# list(lcData.keys())

# %% [markdown]
# (I used `list()` above because Jupyter renders it a bit more nicely than the `dict_keys` object). There is a lot to take in there, but most of those entries are just light curve data.
# 
# Let's first just check the keys that gave me some information about the light curve generically:

# %%
# print(f"Binning: {lcData['Binning']}")
# print(f"TimeFormat: {lcData['TimeFormat']}")
# print(f"T0: {lcData['T0']}")

# %% [markdown]
# As a quick aside, you may wonder why we bother with this array, but we can't just step over the keys in `lcData` - 'Binning', for example, is not a light curve, and maybe in the future we'll want to add other things, so 'Datasets' is handy. 
# 
# There are a lot of datasets in this example, because we set both `incbad` and `nosys` to "both", so we got all data with/out missing centroids and with/out WT-mode systematics (if you don't know what I'm talking about, see the [the light curve documentation](https://www.swift.ac.uk/user_objects/lc_docs.php#systematics).
# 
# 
# The contents of the datasets were discussed in [light curve `dict` documentation](https://www.swift.ac.uk/API/ukssdc/structures.md#the-light-curve-dict) (I'm sounding like a broken record, I know), so I'm not going to spend time on it here, except to show you one entry as an example:

# %%
# lcData['PC']

# %% [markdown]
# (If you're reading the markdown version, I have probably truncated the output, removed some of the rows).
# 
# OK, so that's what `returnData=True` does. If we supply a list of GRBs then, as you should expect, we get a `dict` of light curves `dict`s; the top level is indexed either by GRB name, or targetID, depending on how you called the function, so:
# 

# %%
# lcData = udg.getLightCurves(GRBName=("GRB 211211A","GRB 230328B"),
#                             saveData=False,
#                             returnData=True)

# %%
# list(lcData.keys())

# %% [markdown]
# I trust this doesn't come as a surprise! Nor should the fact that each of these in turn is a light curve `dict` similar to that above (although this time I left `nosys` and `incbad` to their defaults). I can prove this easily enough:

# %%
# list(lcData['GRB 211211A'].keys())

# %% [markdown]
# #### Plotting light curves
# 
# If we've downloaded a light curve then we can make use of the [module-level `plotLightCurve()` function](https://www.swift.ac.uk/API/ukssdc/commonFunc.md#plotlightcurve) to give us a quick plot. I'm not going to repeat the `plotLightCurve()` documentation here, but I will note that its first argument is a single light curve `dict` so if, as in our case here, we downloaded multiple GRB light curves, we have to provide one of them to the function, like this:

# %%
# import swifttools.ukssdc.data.GRB as udg
# from swifttools.ukssdc import plotLightCurve
# lcData = udg.getLightCurves(GRBName=("GRB 211211A","GRB 230328B"),
#                             saveData=False,
#                             returnData=True)
# fig, ax = plotLightCurve(lcData['GRB 211211A'],
#                          whichCurves=('WT_incbad', 'PC_incbad'),
#                          xlog=True,
#                          ylog=True
#                         )

# %% [markdown]
# You'll note I captured the return in variable names which will be familiar to users of `pyplot` and can be ignored by everyone else because you'll need to be familiar with `pyplot` to take advantage of the fact that `plotLightCurve` returnes them for you.

# %% [markdown]
# <a id='ban'></a>
# ## Burst analyser
# 
# This API gives access to all the data in the burst analyser, with some enhancements to which I will return in a moment. First, a reminder that this webpage is documenting the API, not the burst analyser, so if you don't understand some of what I'm discussing, I advise you to look at the [burst analyser paper](https://ui.adsabs.harvard.edu/abs/2010A%26A...519A.102E/abstract) and/or [online documentation](https://www.swift.ac.uk/burst_analyser/docs.php). Surprisingly enough, we get at burst analyser data with the function: `getBurstAnalyser()`. The burst analyser is built on top of light curves, and I will in places refer to things in the [light curve section](#light-curves) so it may be advisable to read that before this.
# 
# The burst analyser is a rather complex beast, because it has so many different datasets in it. We have 3 instruments, multiple energy bands, unabsorbed and observed fluxes, hardness ratios, inferred photon indices and energy conversion factors... oh yes, and a few different options for the BAT binning and multiple UVOT filters. All in all, it's complicated. On the website this is all managed through dividing the page into sections and giving various controls. For the API, it's handled by defining a data structure, the burst analyser `dict`, that contains everything, allowing you to explore it. This does mean that there are an awful lot of parameters available for us to consider when saving burst analyser data.
# 
# We will return to the burst analyser `dict`, and those parameters, in a minute, but first let's discuss the concept of saving data, because for the burst analyser this is a little more complicated than for the above products.
# 
# Conceptually, the situation is exactly the same as for light curves and spectra: you can either get the files from the website and save them straight to disk, or you can download them into a `dict` and, if you want to, write files to disk based on that. However, the way the files are organised for the website is optimised for online visualisation rather than access and manipulation, whereas the whole point of the API *is* access and manipulation. As a result, I decided to prioritise usefulness over uniformity, and give `getBurstAnalyser()` *three* options for what it does:
# 
# * `returnData` - (default `False`), whether or not the function should return the burst analyser `dict`.
# * `saveData` - (default: `True`), whether the data in the burst analyser `dict` to be saved to disk.
# * `downloadTar` - (default: `True`), whether the burst analyser tar file should be saved to disk.
# 
# Here `saveData` is effectively "the third way" defined for light curves and spectra, but automated: the data are downloaded into a burst analyser `dict` and then saved from that (that `dict` is, however, discarded unless `returnData=True`), but the files saved are, I think, much more helpful than those you would get just by grabbing the files from the website. `downloadTar` does let you pull the files straight from the web, and has accompanying boolean arguments `extract` and `removeTar` which lets you, er, extract the data from the `tar` file and then remove said file. I haven't included an explicit demonstration of `downloadTar=True` here because it should be obvious what I mean, it's easy for you to test, and it creates a *lot* of files; *and* because I personally advocate `saveData=True` instead. Oh, and as with the other products, these three parameters are not mutually exclusive, then can all be `True` (or `False` bur I still can't see why you would do that).
# 
# Before we plunge into some demos, though, I should elaborate briefly on the above: what are the 'enhancements' I referred to, and the difference between the files that you get with `downloadTar` and `saveData`? There are two parts to this.
# 
# First: most of the light curves in the burst analyser data actually consist of three time series: the flux (i.e. the light curve), the photon index and the energy conversion factor (ECF), and on the website (and in the downloable `tar` file) these are all in separate files, even though they share a common time axis. So if you want to do any manipulation, you have to read multiple files and then join them on said time axis. With `saveData=True`, this is done for you, so for each light curve you get **one** file that has columns for flux, photon index, ECF etc.
# 
# Second: The issue of error propagation for the burst analyser is complicated (do see [burst analyser paper](https://ui.adsabs.harvard.edu/abs/2010A%26A...519A.102E/abstract) and/or [online documentation](https://www.swift.ac.uk/burst_analyser/docs.php) if you want details). As the documentation explains, in the light curves online (and in the `tar` file), the errors on the flux values are derived solely from the error on the underlying count-rate light curves, the uncertainty in the spectral shape and hence flux conversion are **not** propagaged into those errors. The reasons for this are subtle (but important, and discussed in the documentation), and of course you can do this propagation yourself if you grab all the files. However, in the API, we do this for you! The downloaded data (which we will explore soon) contain two sets of errors &#8212; with and without propagation &#8212; and you can choose which to save to disk as we shall demonstrate in a second.
# 
# Right, that's enough talk, let's get to work.

# %% [markdown]
# <a id='ban_dict'></a>
# ### Getting the burst analyser data into a variable
# 
# I'm going to begin the burst analyser tutorial with the `returnData=True` case (unlike for the earlier products) because this introduces the data that we save to disk with `saveData=True`.
# 
# As with all previous products, this needs the `GRBName` or `targetID` arguments (see [the introduction](#intro)) which can be single values or lists/tuples, and it returns a `dict`. Unlike the light curve and spectral `dict`s, the burst analyser `dict` only appears for GRBs, and so while it is described in [the data structure documentation](https://www.swift.ac.uk/API/ukssdc/structures.md#the-burst-analyser-dict) it is only touched on lightly and I will give a full demonstration here.
# 
# The burst analyser `dict` is not too complicated, and in concept is intentionally reminscent of the way the spectral and light curve `dict`s were built. The burst analyser `dict` has (up to three) layers:
# 
# * Instruments
#  * [BAT binning & HR Data]
#    * Light curves & HR Data
#    
# For obvious reasons, the middle layer is only present for the BAT data, and the different instruments have slighly different contents as we'll see. 
# 
# You can see a detailed schematic in [the data structure documentation](https://www.swift.ac.uk/API/ukssdc/structures.md#the-burst-analyser-dict) but let's instead here explore interatively. First, let's get a single GRB:

# %%
# data = udg.getBurstAnalyser(GRBName="GRB 220715B",
#                                 returnData=True,
#                                 # saveData=False

# %% [markdown]
# Right, now we can explore `data`. The top level of this `dict` is all about the instruments:

# %%
# data.keys()

# %% [markdown]
# If you've followed any of the other data structures you can probably guess what this means. The 'Instruments' entry is a list, telling us what intstruments' data we have; the other entries are all the `dict`s containing those data, obviously indexed by the instrument, so:

# %%
# data['Instruments']

# %% [markdown]
# should not be a surprise.
# 
# You will note that, as for the website, the BAT data, and the BAT data without spectral evolution are separate. I spent a while looking into putting them both inside the same entry and then decided it was much more sensible to keep them separate. We'll explore these data, one instrument at a time. The details of this `dict` differ slightly for each instrument so we'll go through them separately:

# %% [markdown]
# #### XRT data
# 
# If you refer [way back up this notebook to the burst analyser `dict` introduction](#ban_dict), you'll remember that the 'binning' layer is only present for BAT, because this is the only instrument for which the burst analyser has multiple binning options. So, when we explore XRT data we should come straight into a light curve `dict`.

# %%
# data = udg.getBurstAnalyser(GRBName="GRB 220715B",
#                                 returnData=True,
#                                 saveData=False,
#                            instruments=('XRT', 'UVOT'))
# list(data['XRT'].keys())

# %% [markdown]
# And we do! Although eagle-eyed readers will realise there is an extra "HRData_PC" entry. This is analogous to the BAT entry, giving the hardness ratio and its conversion to photon index and ECF. The HR data is separated out for the two XRT modes; although this GRB only has PC mode data. For completeness, let's have a quick look at things:

# %%
# data['XRT']['HRData_PC']

# %% [markdown]
# I must point out one little problem here: there should, really, be an 'HRData_incbad' key here, and there isn't, which is why there is only one bin in this hardness ratio. It turns out that this file doesn't exist (the \_incbad hardness ratio exists *and is used to create the burst analyser light curves*, but this nice, combined file of everything isn't saved to disk). I will look into fixing this, and I guess I'll have to remake all the burst analyser data to create all the files(!) which may take a while, but the problem is in the burst analyser, not the API. When I've fixed it, the '\_incbad' entries will appear.
# 
# The other things are just the flux light curves, analogous to the BAT one:

# %%
# print(data['XRT']['Density_PC_incbad'].to_string())

# %%
# print(data['XRT']['Density_WT_incbad'].to_string())

# %%
# print(data['XRT']['ObservedFlux_PC_incbad'].to_string())

# %%
# print(data['XRT']['XRTBand_PC_incbad'].to_string())

# %%
# # print(data['XRT']['Density_WT_incbad']['Time'][28])
# # print(data['XRT']['Density_PC_incbad']['Time'][0])
# # print(data['XRT']['Density_WT_incbad']['Flux'][28])
# # print(data['XRT']['Density_PC_incbad']['Flux'][0])
# print(data['XRT']['Density_PC_incbad']['Time'][0])
# # print(data['XRT']['Density_WT_incbad']['Time'][0])

# %%
# import swifttools.ukssdc.data.GRB as udg
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# data = udg.getBurstAnalyser(GRBName="GRB 080727C",
#                                 returnData=True,
#                                 saveData=False,
#                            instruments=('XRT', 'UVOT'))
# # print(list(data['XRT'].keys()))
# # Flux_dens_pos=pd.DataFrame(np.concatenate((data['XRT']['Density_WT_incbad']['FluxPos'],\
# #                         data['XRT']['Density_PC_incbad']['FluxPos'])))
# # Flux_dens_neg=pd.DataFrame(np.concatenate((data['XRT']['Density_WT_incbad']['FluxNeg'],\
# #                         data['XRT']['Density_PC_incbad']['FluxNeg'])))
# # time_pos=pd.DataFrame(np.concatenate((data['XRT']['Density_WT_incbad']['TimePos'],\
# #                         data['XRT']['Density_PC_incbad']['TimePos'])))
# # time_neg=pd.DataFrame(np.concatenate((data['XRT']['Density_WT_incbad']['TimeNeg'],\
# #                         data['XRT']['Density_PC_incbad']['TimeNeg'])))
# # time=np.concatenate((data['XRT']['Density_WT_incbad']['Time'], \
# #                 data['XRT']['Density_PC_incbad']['Time']))
# # flux=np.concatenate((data['XRT']['Density_WT_incbad']['Flux'], \
# #                 data['XRT']['Density_PC_incbad']['Flux']))
# # flux_error=np.transpose(np.column_stack((Flux_dens_pos.fillna(0),\
# #         np.abs(Flux_dens_neg.fillna(0)))))
# # time_error=np.transpose(np.column_stack((time_pos.fillna(0),\
# #         np.abs(time_neg.fillna(0)))))
# Flux_dens_pos=data['XRT']['Density_PC_incbad']['FluxPos']
# Flux_dens_neg=data['XRT']['Density_PC_incbad']['FluxNeg']
# time_pos=data['XRT']['Density_PC_incbad']['TimePos']
# time_neg=data['XRT']['Density_PC_incbad']['TimeNeg']
# time=data['XRT']['Density_PC_incbad']['Time']
# flux=data['XRT']['Density_PC_incbad']['Flux']
# flux_error=np.transpose(np.column_stack((Flux_dens_pos.fillna(0),\
#         np.abs(Flux_dens_neg.fillna(0)))))
# time_error=np.transpose(np.column_stack((time_pos.fillna(0),\
#         np.abs(time_neg.fillna(0)))))
# # Flux_dens_pos=data['XRT']['Density_WT_incbad']['FluxPos']
# # Flux_dens_neg=data['XRT']['Density_WT_incbad']['FluxNeg']
# # time_pos=data['XRT']['Density_WT_incbad']['TimePos']
# # time_neg=data['XRT']['Density_WT_incbad']['TimeNeg']
# # time=data['XRT']['Density_WT_incbad']['Time']
# # flux=data['XRT']['Density_WT_incbad']['Flux']
# # flux_error=np.transpose(np.column_stack((Flux_dens_pos.fillna(0),\
# #         np.abs(Flux_dens_neg.fillna(0)))))
# # time_error=np.transpose(np.column_stack((time_pos.fillna(0),\
# #         np.abs(time_neg.fillna(0)))))
# # z=2.3739
# # timeshift=200*(1+z)
# # # timeshift=2000
# # print(timeshift)
# # start_time=np.where(time >= timeshift)[0][0]
# # print(time[start_time])
# timeshift=10**(-0.301)*36600
# # timeshift=0
# start_time=(np.where(time >= timeshift)[0])[0]
# timeclaw=10**(0.301)*36600
# stop_time=max(np.where(time<= timeclaw)[0])
# coefficients = np.polyfit(np.log10(time[start_time:stop_time]), \
#                           np.log10(flux[start_time:stop_time]), 1)
# polynomial = np.poly1d(coefficients)
# log10_y_fit = polynomial(np.log10(time[start_time:stop_time]))
# plt.errorbar(x=time, \
#              y=flux, \
#             xerr=time_error, yerr=flux_error, ls='none')
# plt.plot(time[start_time:stop_time], 10**log10_y_fit, '--')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'Time (s)', \
#            color='0.8', size=14)
# plt.ylabel(r'${Flux}_{XRT}$ (0.3-10 keV) $(erg*cm^{-2}*s^{-1})$', \
#            color='0.8', size=14)
# plt.tick_params(axis='x', colors='0.8')
# plt.tick_params(axis='y', colors='0.8')
# # print(polynomial[1])
# print(time[start_time:stop_time])


