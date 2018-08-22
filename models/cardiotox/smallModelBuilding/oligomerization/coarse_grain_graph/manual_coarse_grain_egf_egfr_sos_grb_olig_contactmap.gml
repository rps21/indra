Creator	"yFiles"
Version	"2.14"
graph
[
	hierarchic	1
	label	""
	directed	1
	node
	[
		id	0
		label	"EGF"
		graphics
		[
			x	35.0
			y	92.0
			w	30.0
			h	30.0
			type	"roundrectangle"
			fill	"#D2D2D2"
			outline	"#000000"
			topBorderInset	0.0
			bottomBorderInset	0.0
			leftBorderInset	0.0
			rightBorderInset	0.0
		]
		LabelGraphics
		[
			text	"EGF"
			fontSize	12
			fontStyle	"bold"
			fontName	"Dialog"
			anchor	"t"
		]
		isGroup	1
	]
	node
	[
		id	1
		label	"EGFR"
		graphics
		[
			x	141.12222222222223
			y	92.0
			w	30.0
			h	30.0
			type	"roundrectangle"
			fill	"#D2D2D2"
			outline	"#000000"
			topBorderInset	0.0
			bottomBorderInset	0.0
			leftBorderInset	0.0
			rightBorderInset	0.0
		]
		LabelGraphics
		[
			text	"EGFR"
			fontSize	12
			fontStyle	"bold"
			fontName	"Dialog"
			anchor	"t"
		]
		isGroup	1
	]
	node
	[
		id	2
		label	"GRB2"
		graphics
		[
			x	231.1579365079365
			y	92.0
			w	30.0
			h	30.0
			type	"roundrectangle"
			fill	"#D2D2D2"
			outline	"#000000"
			topBorderInset	0.0
			bottomBorderInset	0.0
			leftBorderInset	0.0
			rightBorderInset	0.0
		]
		LabelGraphics
		[
			text	"GRB2"
			fontSize	12
			fontStyle	"bold"
			fontName	"Dialog"
			anchor	"t"
		]
		isGroup	1
	]
	node
	[
		id	3
		label	"SOS1"
		graphics
		[
			x	189.89007936507937
			y	15.0
			w	30.0
			h	30.0
			type	"roundrectangle"
			fill	"#D2D2D2"
			outline	"#000000"
			topBorderInset	0.0
			bottomBorderInset	0.0
			leftBorderInset	0.0
			rightBorderInset	0.0
		]
		LabelGraphics
		[
			text	"SOS1"
			fontSize	12
			fontStyle	"bold"
			fontName	"Dialog"
			anchor	"t"
		]
		isGroup	1
	]
	edge
	[
		source	0
		target	1
		graphics
		[
			fill	"#000000"
			sourceArrow	"delta"
			Line
			[
				point
				[
					x	35.0
					y	92.0
				]
				point
				[
					x	80.12222222222222
					y	92.0
				]
				point
				[
					x	80.12222222222222
					y	99.5
				]
				point
				[
					x	141.12222222222223
					y	92.0
				]
			]
		]
		edgeAnchor
		[
			xSource	1.0
			xTarget	-1.0
			yTarget	0.5
		]
	]
	edge
	[
		source	1
		target	1
		graphics
		[
			fill	"#000000"
			sourceArrow	"delta"
			Line
			[
				point
				[
					x	141.12222222222223
					y	92.0
				]
				point
				[
					x	110.62222222222223
					y	84.5
				]
				point
				[
					x	110.62222222222223
					y	61.5
				]
				point
				[
					x	133.62222222222223
					y	61.5
				]
				point
				[
					x	141.12222222222223
					y	92.0
				]
			]
		]
		edgeAnchor
		[
			xSource	-1.0
			ySource	-0.5
			xTarget	-0.5
			yTarget	-1.0
		]
	]
	edge
	[
		source	1
		target	3
		graphics
		[
			fill	"#000000"
			sourceArrow	"delta"
			Line
			[
				point
				[
					x	141.12222222222223
					y	92.0
				]
				point
				[
					x	148.62222222222223
					y	45.5
				]
				point
				[
					x	182.39007936507937
					y	45.5
				]
				point
				[
					x	189.89007936507937
					y	15.0
				]
			]
		]
		edgeAnchor
		[
			xSource	0.5
			ySource	-1.0
			xTarget	-0.5
			yTarget	1.0
		]
	]
	edge
	[
		source	1
		target	2
		graphics
		[
			fill	"#000000"
			sourceArrow	"delta"
		]
		edgeAnchor
		[
			xSource	1.0
			xTarget	-1.0
		]
	]
	edge
	[
		source	2
		target	3
		graphics
		[
			fill	"#000000"
			sourceArrow	"delta"
			Line
			[
				point
				[
					x	231.1579365079365
					y	92.0
				]
				point
				[
					x	231.1579365079365
					y	45.5
				]
				point
				[
					x	197.39007936507937
					y	45.5
				]
				point
				[
					x	189.89007936507937
					y	15.0
				]
			]
		]
		edgeAnchor
		[
			ySource	-1.0
			xTarget	0.5
			yTarget	1.0
		]
	]
]
