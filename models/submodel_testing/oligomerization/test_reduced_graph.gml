graph [
node [
id 0 label "EGF" isGroup 1  graphics [ type "roundrectangle" fill "#D2D2D2" outline "#000000"  ] LabelGraphics [ text "EGF" anchor "t" fontStyle "bold"  ]
]
node [
id 3 label "EGFR" isGroup 1  graphics [ type "roundrectangle" fill "#D2D2D2" outline "#000000"  ] LabelGraphics [ text "EGFR" anchor "t" fontStyle "bold"  ]
]
node [
id 14 label "GRB2" isGroup 1  graphics [ type "roundrectangle" fill "#D2D2D2" outline "#000000"  ] LabelGraphics [ text "GRB2" anchor "t" fontStyle "bold"  ]
]
node [
id 24 label "SOS1" isGroup 1  graphics [ type "roundrectangle" fill "#D2D2D2" outline "#000000"  ] LabelGraphics [ text "SOS1" anchor "t" fontStyle "bold"  ]
]
edge [ source 3 target 3  graphics [ fill "#000000"  ] ]
edge [ source 14 target 3  graphics [ fill "#000000"  ] ]
edge [ source 3 target 24  graphics [ fill "#000000"  ] ]
edge [ source 0 target 3  graphics [ fill "#000000"  ] ]
edge [ source 14 target 0  graphics [ fill "#000000"  ] ]
edge [ source 14 target 24  graphics [ fill "#000000"  ] ]
]
