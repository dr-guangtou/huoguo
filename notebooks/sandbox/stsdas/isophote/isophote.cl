procedure isophote()
string	mode="al"

begin
	package isophote

	task ellipse = "isophote$x_isophote.e"
	task model   = "isophote$x_isophote.e"
	task map     = "isophote$x_isophote.e"
	task isoexam = "isophote$x_isophote.e"
	task isoplot = "isophote$isoplot.cl"
	task isopall = "isophote$isopall.cl"
	task isomap  = "isophote$isomap.cl"
	task isoimap = "isophote$isoimap.cl"
	task bmodel  = "isophote$bmodel.cl"
	hidetask	model
	hidetask	map

	# Psets
	task geompar     = "isophote$geompar.par"
	task controlpar  = "isophote$controlpar.par"
	task samplepar   = "isophote$samplepar.par"
	task magpar      = "isophote$magpar.par"

	cl()
end
