from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JetData-RunB2017-v1'
config.General.workArea = 'JetData-RunB2017-v1'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName  = 'flatData-nanoNtuple-cfg.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2017B-31Mar2018-v1/MINIAOD'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.runRange = '294927-306462'
config.Data.outputDatasetTag = 'CRAB3_JetData-RunB2017'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"
