from PhysicsTools.TensorFlow.tfGraphDefProducer_cfi import tfGraphDefProducer as _tfGraphDefProducer
superclusteringTf = _tfGraphDefProducer.clone(
    ComponentName = "superclusteringTf",
    FileName = "RecoHGCal/TICL/data/tf_models/supercls_v1.pb"
)
