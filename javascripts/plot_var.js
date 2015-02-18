/**
 * @author Yafeng Hao, Liran Juan
 */
var colD = "#d73027"; //denovo mut color
var colD_hover = "#f46d43"; //denovo mut color hover
var colP = "#4575b4"; //paternal var color
var colP_hover = "#74add1"; //paternal var color hover
var colM = "#1a9850"; //maternal var color
var colM_hover = "#66bd63"; //maternal var color
var colO = "#000000"; //black
var colO_hover = "#636363"; //black hover
var colO_hover2 = "#bdbdbd"; //black hover
var colInter = "#542788";//intersection color
var colInter_hover = "#998ec3";//intersection color hover
var colDiff = "#542788";//difference color
var colDiff_hover = "#998ec3";//difference color hover

///////////////
var restrictX;
var restrictY;
var gdf_list;
var gdf_icon;

var lds_flag = 0;
var lds_btnset;
var parents_area;
var csi_text;
var compared_indiv_area;

var legendset;
var legendset_share_diff;
var blueobj,lg_blue,blackobj,lg_black;

var ncbi_url = "http://www.ncbi.nlm.nih.gov/nuccore/";
var hgnc_url = "http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:";

var bin_set;

var warning_set;
var bins = [];
var rightORleft;

var vdetail = [];
var floatboxset;
var floatbox = [];

var recom_GB = [];
var recom_set;
///////////////

var functional_vnum = 0;
var all_lds ;
function init_genes_omim(){
	variants = [];
	variants_byid = {};
	functional_v = {};
	functional_vPointer = [];
	genes = [];
	symbols = [];
	css = 0;
	cst = 0;
	cssObj = null;
	cstObj = null;

	all_lds = R.set();

	symbols[0] = "ALL";
	functional_vnum = 0;
	var symbol_map = {}
	var sPointer = 1;
	

	var esNodes = req3.responseXML.getElementsByTagName(xmlTagElements);
	var enode_start = 0;

	

	for(var i = 0 ; i < esNodes.length ; i++){
		

		var eNodes = esNodes[i].getElementsByTagName(xmlTagElement);
		for(var j = 0 ; j < eNodes.length ; j++){
			if(i == 0){
				genes[j] = {};
				genes[j].id = eNodes[j].getAttribute(xmlAttributeId);
				genes[j].symbol = eNodes[j].getAttribute("Symbol");
				genes[j].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				genes[j].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				genes[j].strand = eNodes[j].getElementsByTagName(xmlTagDirection)[0].firstChild.nodeValue;
				var subNodes = eNodes[j].getElementsByTagName(xmlTagSubElement);
				if(subNodes.length > 0){
					genes[j].subs = [];
					for(var k = 0 ; k < subNodes.length ; k++){
						genes[j].subs[k] = {};
						genes[j].subs[k].type = subNodes[k].getAttribute(xmlAttributeType);
						genes[j].subs[k].from = subNodes[k].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
						genes[j].subs[k].to = subNodes[k].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
					}
				}

				var tr_len = 1;
				if(symbol_map[genes[j].symbol] == undefined){
					symbol_map[genes[j].symbol] = sPointer;
					symbols[sPointer] = [];
					symbols[sPointer][0] = genes[j].symbol;
					sPointer++;
				} else {
					tr_len = symbols[symbol_map[genes[j].symbol]].length;
				}
				symbols[symbol_map[genes[j].symbol]][tr_len] = j;
			}
			else if(esNodes[i].getAttribute(xmlAttributeId)=="_OMIM"){
				gdf_list[j] = {};
				gdf_list[j].symbol= eNodes[j].getAttribute("Symbol");
				gdf_list[j].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				gdf_list[j].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				gdf_list[j].id = eNodes[j].getAttribute(xmlAttributeId);
				gdf_list[j].source = eNodes[j].getElementsByTagName("Source")[0].firstChild.nodeValue;
				/*var sPointerrr = symbol_map[eNodes[j].getAttribute("Symbol")];
				var slength = symbols[sPointerrr].length;
				var newindex = j + enode_start;
				symbols[sPointerrr][slength] = newindex;
				genes[j+enode_start] = {};
				genes[j+enode_start].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				genes[j+enode_start].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				genes[j+enode_start].id = eNodes[j].getAttribute(xmlAttributeId);
				genes[j+enode_start].omim = 1;*/
			}
			
		}
	}

}
function init_individual_vars(){
	variants = [];
	variants_byid = {};
	functional_v = {};
	functional_vPointer = [];
	genes = [];
	symbols = [];
	css = 0;
	cst = 0;
	cssObj = null;
	cstObj = null;

	all_lds = R.set();

	symbols[0] = "ALL";
	functional_vnum = 0;
	var symbol_map = {}
	var sPointer = 1;
	var vsNode = req3.responseXML.getElementsByTagName(xmlTagVariants)[0];
	var vNodes = vsNode.getElementsByTagName(xmlTagVariant);
	for(var i = 0 ; i < vNodes.length ; i++){
		variants[i] = {};
		variants[i].chr = current_chr;
		variants[i].id = vNodes[i].getAttribute(xmlAttributeId);
		variants[i].dd = vNodes[i].getAttribute("dd");

		variants[i].genotype = vNodes[i].getAttribute(xmlAttributeGenotype);
		var splitter = "|";
		if(variants[i].genotype.indexOf("/")>0){
			splitter = "/";
		}
		variants[i].genotypes = variants[i].genotype.split(splitter);

		if(i > 0 && variants[i-1].genotypes[0]*variants[i-1].genotypes[1] > 1 
		&& variants[i].genotypes[0]*variants[i].genotypes[1] > 1
		&& (variants[i].dd == undefined && variants[i-1].dd == undefined
		|| variants[i].dd != undefined && variants[i-1].dd != undefined 
		&& variants[i].dd == variants[i-1].dd)){

			variants[i-1].genotype = variants[i-1].genotype + "(" + variants[i-1].genotypes[0] + ")";
			variants[i].genotype = variants[i].genotype + "(" + variants[i].genotypes[1] + ")";
			variants[i-1].genotypes[1] = 0;
			variants[i].genotypes[0] = 0;
		}
		variants[i].type = vNodes[i].getAttribute(xmlAttributeType);

		if(variants[i].id == undefined || variants[i].id == "." || variants[i].id ==""){
			variants[i].id = "Unnamed "+variants[i].type;
		}
		variants[i].from = parseInt(vNodes[i].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue);
		variants[i].to = parseInt(vNodes[i].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue);
		variants[i].selected = false;

		if(vNodes[i].getElementsByTagName(xmlTagLetter).length > 0){
			variants[i].letter = vNodes[i].getElementsByTagName(xmlTagLetter)[0].childNodes[0].nodeValue;
		} else {
			variants[i].letter = "--";
		}

		variants_byid[current_chr+":"+variants[i].from+":"+variants[i].to+":"+(variants[i].letter=="--"?"-":variants[i].letter)] = i;
		variants[i].paternal = "Y";
		variants[i].maternal = "Y";
		variants[i].functional = "--";
	}

	var esNodes = req3.responseXML.getElementsByTagName(xmlTagElements);
	var enode_start = 0;

	function match_functional_variants(vPointer, from, to, letter, id){
		var letter_trans = letter;

		if(letter == "(" || letter == ")"){
			letter_trans = "ASS/DSS loss";
		}else if (letter == "#"){
			letter_trans = "Frame shifting";
		}else if (letter == "^"){
			letter_trans = "Initiator loss";
		}else{
			var texts = letter.split(":");
			if(texts[0].indexOf("$") >= 0 && texts[1].indexOf("$") >= 0){
				letter_trans = "Stop -> Stop";
			}else if(texts[0].indexOf("$") >= 0){
				letter_trans = "Stop loss";
			}else if(texts[1].indexOf("$") >= 0){
				letter_trans = "Stop gain";
			}else {
				letter_trans = texts[0] + "->" + texts[1];
			}
		}

		while(vPointer > 0 && variants[vPointer].to > to){
			vPointer--;
		}
		while(vPointer < variants.length-1 && variants[vPointer].from < from){
			vPointer++;
		}
		if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
			if(functional_v[from+":"+to+":"+letter] == undefined){
				functional_v[from+":"+to+":"+letter] = vPointer;
				functional_vnum++;
				var fvidx = 0; 
				if(variants[vPointer].functional == "--"){
					//variants[vPointer] has function for the first time
					variants[vPointer].functional = letter_trans;
					variants[vPointer].functional_v = [];
					functional_vPointer[functional_vPointer.length] = vPointer;
				}else{
					//variants[vPointer] may play different roles in different allele/transcript
					variants[vPointer].functional += ","+letter_trans;
					fvidx = variants[vPointer].functional_v.length;
				}
				variants[vPointer].functional_v[fvidx] = {};
				variants[vPointer].functional_v[fvidx].from = from;
				variants[vPointer].functional_v[fvidx].to = to;
				variants[vPointer].functional_v[fvidx].letter = letter;
			}
		}else if(from <= variants[vPointer].from && to >= variants[vPointer].to){
			while(vPointer > 0 && from <= variants[vPointer-1].from && to >= variants[vPointer-1].to){
				vPointer--;
			}
			while(vPointer < variants.length && from <= variants[vPointer].from && to >= variants[vPointer].to){
				if(id == variants[vPointer].id){
					if(functional_v[from+":"+to+":"+letter] == undefined){
						functional_v[from+":"+to+":"+letter] = vPointer;
						functional_vnum++;
						var fvidx = 0; 
						if(variants[vPointer].functional == "--"){
							//variants[vPointer] has function for the first time
							variants[vPointer].functional = letter_trans;
							variants[vPointer].functional_v = [];
							functional_vPointer[functional_vPointer.length] = vPointer;
						}else{
							//variants[vPointer] may play different roles in different allele/transcript
							variants[vPointer].functional += ","+letter_trans;
							fvidx = variants[vPointer].functional_v.length;
						}
						variants[vPointer].functional_v[fvidx] = {};
						variants[vPointer].functional_v[fvidx].from = from;
						variants[vPointer].functional_v[fvidx].to = to;
						variants[vPointer].functional_v[fvidx].letter = letter;
					}
					return vPointer;
				} else {
					vPointer++;
				}
			}
			return vPointer-1;
		}
		return vPointer;
	}

	for(var i = 0 ; i < esNodes.length ; i++){
		var veNodes = esNodes[i].getElementsByTagName(xmlTagVariant);
		var vPointer = 0;
		for(var j = 0 ; j < veNodes.length ; j++){
			var from = veNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
			var to = veNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
			var letter = veNodes[j].getElementsByTagName(xmlTagLetter)[0].childNodes[0].nodeValue;
			var id = veNodes[j].getAttribute(xmlAttributeId);
			var type = veNodes[j].getAttribute(xmlAttributeType);

			if(id == "."){
				id = "Unnamed "+type;
			}
			if(from.indexOf(";")<0 && to.indexOf(";")<0){
				match_functional_variants(vPointer, from, to, letter, id);
			} else {
				var froms = from.split(";");
				var tos = to.split(";");
				if(type == "DEL" && tos[1] - froms[0] < 50){
					match_functional_variants(vPointer, froms[0], tos[1], letter, id);
				}else{
					for(var junc_var = 0 ; junc_var < froms.length ; junc_var++){
						if(froms[junc_var]>=current_start && tos[junc_var]<=current_end){
							vPointer = match_functional_variants(vPointer, froms[junc_var], tos[junc_var], letter, id);
						}
					}
				}
			}
		}

		var eNodes = esNodes[i].getElementsByTagName(xmlTagElement);
		for(var j = 0 ; j < eNodes.length ; j++){
			if(i == 0){
				genes[j] = {};
				genes[j].id = eNodes[j].getAttribute(xmlAttributeId);
				genes[j].symbol = eNodes[j].getAttribute("Symbol");
				genes[j].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				genes[j].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				genes[j].strand = eNodes[j].getElementsByTagName(xmlTagDirection)[0].firstChild.nodeValue;
				var subNodes = eNodes[j].getElementsByTagName(xmlTagSubElement);
				if(subNodes.length > 0){
					genes[j].subs = [];
					for(var k = 0 ; k < subNodes.length ; k++){
						genes[j].subs[k] = {};
						genes[j].subs[k].type = subNodes[k].getAttribute(xmlAttributeType);
						genes[j].subs[k].from = subNodes[k].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
						genes[j].subs[k].to = subNodes[k].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
					}
				}

				var tr_len = 1;
				if(symbol_map[genes[j].symbol] == undefined){
					symbol_map[genes[j].symbol] = sPointer;
					symbols[sPointer] = [];
					symbols[sPointer][0] = genes[j].symbol;
					sPointer++;
				} else {
					tr_len = symbols[symbol_map[genes[j].symbol]].length;
				}
				symbols[symbol_map[genes[j].symbol]][tr_len] = j;
			}
			else if(esNodes[i].getAttribute(xmlAttributeId)=="_OMIM"){
				gdf_list[j] = {};
				gdf_list[j].symbol= eNodes[j].getAttribute("Symbol");
				gdf_list[j].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				gdf_list[j].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				gdf_list[j].id = eNodes[j].getAttribute(xmlAttributeId);
				gdf_list[j].source = eNodes[j].getElementsByTagName("Source")[0].firstChild.nodeValue;
				/*var sPointerrr = symbol_map[eNodes[j].getAttribute("Symbol")];
				var slength = symbols[sPointerrr].length;
				var newindex = j + enode_start;
				symbols[sPointerrr][slength] = newindex;
				genes[j+enode_start] = {};
				genes[j+enode_start].from = eNodes[j].getElementsByTagName(xmlTagFrom)[0].firstChild.nodeValue;
				genes[j+enode_start].to = eNodes[j].getElementsByTagName(xmlTagTo)[0].firstChild.nodeValue;
				genes[j+enode_start].id = eNodes[j].getAttribute(xmlAttributeId);
				genes[j+enode_start].omim = 1;*/
			}
			var vNodes = eNodes[j].getElementsByTagName(xmlTagVariant);
			if(vNodes.length > 0){
				if(genes[j].vars == undefined){
					genes[j].vars = {};
				}
				for(var k = 0 ; k < vNodes.length ; k++){
					var from = vNodes[k].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
					var to = vNodes[k].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
					var letter = vNodes[k].getElementsByTagName(xmlTagLetter)[0].childNodes[0].nodeValue;
					if(from.indexOf(";")>0 || to.indexOf(";")>0){
						var froms = from.split(";");
						var tos = to.split(";");
						if(functional_v[froms[0] + ":" + tos[0] + ":" + letter] != undefined){
							from = froms[0];
							to = tos[0];
						} else if (functional_v[froms[1] + ":" + tos[1] + ":" + letter] != undefined){
							from = froms[1];
							to = tos[1];
						} else if (functional_v[froms[0] + ":" + tos[1] + ":" + letter] != undefined){
							from = froms[0];
							to = tos[1];
						}
					}
					if(variants[functional_v[from + ":" + to + ":" + letter]] != undefined){
						for(var fvidx = 0 ; fvidx < variants[functional_v[from+":"+to+":"+letter]].functional_v.length ; fvidx++){
							if(from == variants[functional_v[from+":"+to+":"+letter]].functional_v[fvidx].from 
							&& to == variants[functional_v[from+":"+to+":"+letter]].functional_v[fvidx].to 
							&& letter == variants[functional_v[from+":"+to+":"+letter]].functional_v[fvidx].letter){
								genes[j].vars[from+":"+to+":"+letter] = fvidx;
							}
						}
					}else{
//						alert(from + ":" + to + ":" + letter);
					}
				}
			}
		}
	}
	functional_vPointer.sort(sortNumber);
}

function show_axis(){
	close_shared_different();
	close_select_method();
	compared_indiv_area = R.set();
	var l = R_height - R_top - R_bottom;
	var w = R_width - R_left - R_right;
	var font_size2_text = "14px \"Trebuchet MS\", Arial, sans-serif";
	var axis_set = R.set();
	var axis = R.path("M"+R_left+","+R_top+" L"+R_left+","+(R_height-R_bottom)).attr({stroke:colO,"stroke-width":1});
	var cali = [];
	var warningtext1 = R.text(0,R_top-40,"The requested scope is too large.").attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colD});
	var warningtext2 = R.text(0,R_top-25,"Please wait ...").attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colD});
	warning_set = R.set();
	warning_set.push(warningtext1, warningtext2);
	if(current_end - current_start < 3000000){
		warning_set.hide();
	}else{
		warning_set.show();
	}
	var t_bot, t_top;
	t_top = R.text(0,R_top-10,current_chr+":"+current_start).attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colO});
	t_bot = R.text(0,R_height-R_bottom+10,current_chr+":"+current_end).attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colO});
	axis_set.push(axis, t_bot, t_top);
	ca = 0;
	axis_set.push(cali[ca]);
	if(current_end - current_start > 424){
		for(var ca = 0 ; ca <= 100 ; ca++){
			if(ca%10 == 0){
				cali[ca] = R.path("M"+(R_left-20)+","+(R_top+ca/100*l)+" L"+R_left+","+(R_top+ca/100*l)).attr({stroke:colO,"stroke-width":1});
			}else{
				cali[ca] = R.path("M"+(R_left-10)+","+(R_top+ca/100*l)+" L"+R_left+","+(R_top+ca/100*l)).attr({stroke:colO,"stroke-width":1});
			}
			axis_set.push(cali[ca]);
		}
		axis_set.push(cali[ca]);
	}else{
		if(seqstr != undefined && seqstr != ""){
			axis.attr({stroke:"#FFFFFF"});
			var unit_height = ((R_height-R_bottom-R_top)/seqstr.length).toFixed(3);
			var seqchar = seqstr.split("");
			var unit_y = R_top;
			for(var i = 0; i < seqchar.length; i++){
				var temp_unitobj = R.rect(R_left-15, unit_y, 25, unit_height, 0);
				if(unit_height > 15){
					R.text(R_left-2, unit_y + unit_height/2, seqchar[i]).attr({font:"8px  \"Trebuchet MS\", Arial, sans-serif"});
				}
				unit_y = parseFloat(unit_y) + parseFloat(unit_height);
				switch(seqchar[i]){
					case "A":
						temp_unitobj.attr({stroke:"rgb(249,194,56)",fill:"rgb(249,194,56)"});
					break;
					case "T":
						temp_unitobj.attr({stroke:"rgb(133,122,185)",fill:"rgb(133,122,185)"});
					break;
					case "G":
						temp_unitobj.attr({stroke:"rgb(122,197,131)",fill:"rgb(122,197,131)"});
					break;
					case "C":
						temp_unitobj.attr({stroke:"rgb(236,95,75)",fill:"rgb(236,95,75)"});
					break;
					case "N":
						temp_unitobj.attr({stroke:"rgb(163,165,167)",fill:"rgb(163,165,167)"});
					break;
				}
			}
		}
	}
	//axis_set.hide();
	
	var brwplotCanvas = $(document.getElementById("brwplot"));
	//var touch_x = restrictX - brwplotCanvas.position().left;
	//var touch_y = restrictY - brwplotCanvas.position().top;
	var markline = R.path("M"+(R_left-35)+","+R_top+" L"+(R_left+R_width)+","+R_top).attr({"stroke-dasharray": "-",stroke:"#000","stroke-width":1});
	var marklable = R.text(R_left+1, R_top-8, "").attr({font:"8px  \"Trebuchet MS\", Arial, sans-serif"});
	
	var start_y;
	var dis_y = 0;
	var searcharea = R.rect(R_left-30, R_top,R_width+35,1).attr({fill:"#ffe4b5", "stroke-width":0, "fill-opacity":0.4});
	searcharea.hide();
	var	start_line;
	var toucharea = R.rect(R_left-25, R_top - 1, 50, R_height - R_bottom - R_top + 1, 0).attr({stroke: "#000", fill: "#000", opacity: 0, cursor: "move"});
	toucharea.drag(dragmove, dragstart, dragend);
	function dragstart(x,y) {
		start_y = y - brwplotCanvas.position().top;
		searcharea.animate({transform: "t0 "+(start_y-R_top)});
		searcharea.show();
		start_line = R.path("M"+(R_left-35)+","+start_y+" L"+(R_left+R_width)+","+start_y).attr({"stroke-dasharray": "-",stroke:"#000","stroke-width":1});
	}
	function dragmove(dx,dy) {
		if(dy>0){
			if(dy + start_y < R_height - R_bottom){
				searcharea.animate({height: dy});
				dis_y = dy;
			}else{
				searcharea.animate({height: (R_height - R_bottom - start_y)});
				dis_y = R_height - R_bottom - start_y;
			}			
		}else if(dy<0){
			if(dy + start_y > R_top){
				searcharea.animate({transform: "t0 "+(start_y-R_top+dy)});
				searcharea.animate({height: (0-dy)});
				dis_y = dy;
			}else{
				searcharea.animate({transform: "t0 0"});
				searcharea.animate({height: (start_y - R_top)});
				dis_y = R_top - start_y;
			}
		}
	}
	function dragend() {
		start_line.remove();
		searcharea.hide();
		var searchLength = current_end - current_start;
		var index1 = Math.round(searchLength * (start_y-R_top)/(R_height-R_bottom-R_top) + current_start);
		var index2 = Math.round(searchLength * (start_y-R_top+dis_y)/(R_height-R_bottom-R_top) + current_start);
		if(Math.abs(index1-index2) > 2){	
			if(dis_y > 0){
				load_family_genome(current_chr, index1, index2);
			}else{
				load_family_genome(current_chr, index2, index1);
			}
		}else{
			if(dis_y > 0){
				load_family_genome(current_chr, index1, index1+3);
			}else{
				load_family_genome(current_chr, index2, index2+3);
			}
		}
	}



	markline.hide();
	marklable.hide();
	toucharea.mouseover(function(){
		//alert("ok "+ restrictX + " " + restrictY + " " + brwplotCanvas.position().left);
		var touch_x = restrictX - brwplotCanvas.position().left;
		var touch_y = restrictY - brwplotCanvas.position().top;
		var searchLength = current_end - current_start + 1;
		//var refSeqIndex = Math.round(searchLength * (touch_y-R_top-(R_height-R_bottom-R_top)/(searchLength * 2))/(R_height-R_bottom-R_top) + current_start);
		var refSeqIndex = Math.floor((touch_y-R_top)/((R_height-R_bottom-R_top)/searchLength) + current_start);

		markline.stop().animate({transform:"t"+(R_left-25)+" "+(touch_y-R_top)},0);
		marklable.stop().animate({transform:"t"+(R_left+1)+" "+(touch_y-R_top)},0);
		marklable.attr({text:refSeqIndex});
		markline.show();
		marklable.show();
	});

	toucharea.mousemove(function(){
		var touch_x = restrictX - brwplotCanvas.position().left;
		var touch_y = restrictY - brwplotCanvas.position().top;
		var searchLength = current_end - current_start + 1;
		//var refSeqIndex = Math.round(searchLength * (touch_y-R_top)/(R_height-R_bottom-R_top) + current_start);
		var refSeqIndex = Math.floor((touch_y-R_top)/((R_height-R_bottom-R_top)/searchLength) + current_start);
		markline.stop().animate({transform:"t"+(R_left-25)+" "+(touch_y-R_top)},0);
		marklable.stop().animate({transform:"t"+(R_left+1)+" "+(touch_y-R_top)},0);
		marklable.attr({text:refSeqIndex});
		markline.show();
		marklable.show();
	});

	toucharea.mouseout(function(){
		markline.hide();
		marklable.hide();
	});

	if(individuals[csi] != undefined){
		var color = colO;
		var color_text = colO;
		var color_text_hover = colO_hover;
		parents_area = R.set();
		if(individuals[individuals[csi].fid] != undefined && individuals[individuals[csi].mid] != undefined  && individuals[individuals[csi].fid].ifs == "true" && individuals[individuals[csi].mid].ifs == "true"){
			color = "#FFF";
			if(individuals[individuals[csi].fid].affected == "1"){
				color = colO;
				color_hover = colO_hover;
			}
			parents_area.push(R.rect(R_left+w/2+20,5,20,20,0).attr({fill:color,stroke:colO,"stroke-width":1}));
	
			color = "#FFF";
			if(individuals[individuals[csi].mid].affected == "1"){
				color = colO;
				color_hover = colO_hover;
			}
			parents_area.push(R.ellipse(R_left+w/2-30,15,10,10).attr({fill:color,stroke:colO,"stroke-width":1}));
			parents_area.push(R.path("M"+(R_left+w/2-20)+",15 L"+(R_left+w/2+20)+",15").attr({stroke:colO,"stroke-width":1}));
			parents_area.push(R.path("M"+(R_left+w/2)+",15 L"+(R_left+w/2)+","+(R_top-25)).attr({stroke:colO,"stroke-width":1}));
	
			parents_area.push(R.text(R_left+w/2+42,15,individuals[csi].fid).attr({font:font_size2_text,fill:colP,"text-anchor":"start"}));
			parents_area.push(R.text(R_left+w/2-42,15,individuals[csi].mid).attr({font:font_size2_text,fill:colM,"text-anchor":"end"}));
			color = colD;
			color_text = colD;
			color_text_hover = colD_hover;
			
			if(compare_method != "trioAnalysis"){
				parents_area.hide();
			}
		}
		csi_text = R.text(R_left+w/2+15,R_top-15,csi).attr({font:font_size2_text,fill:color,"text-anchor":"start"});

		color = "#FFF";
		if(individuals[csi].affected == "1"){
			color = colO;
			color_hover = colO_hover;
		}
		var indObj = null;
		if(individuals[csi].sex == "1"){
			indObj = R.rect(R_left+w/2-10,R_top-25,20,20,0).attr({fill:color,stroke:colO,"stroke-width":1});
		} else if(individuals[csi].sex == "2"){
			indObj = R.ellipse(R_left+w/2,R_top-15,10,10).attr({fill:color,stroke:colO,"stroke-width":1});
		}else{
			indObj = R.path("M"+(R_left+w/2)
					+" "+(R_top-29)
					+"L"+(R_left+w/2+14)
					+" "+(R_top-15)
					+"L"+(R_left+w/2)
					+" "+(R_top-1)
					+"L"+(R_left+w/2-14)
					+" "+(R_top-15)
					+"Z").attr({fill:color,stroke:colO,"stroke-width":1.5});
		}
		(function(indObj){
			indObj[0].style.cursor = "pointer";
			indObj[0].onmouseover = function(){
				indObj.animate({fill:"#999"},200);
				csi_text.animate({fill:color_text_hover});
				R.safari();
			};
			indObj[0].onmouseout = function(){
				indObj.animate({fill:color},200);
				csi_text.animate({fill:color_text});
				R.safari();
			};
			indObj[0].onclick = function(){
				if(document.getElementById("SD_window").style.display == "none" && 
					document.getElementById("select_method").style.display == "none"){
					if(individuals[csi].fid != "0" && individuals[csi].mid != "0"){
						call_select_method();
					}else{
						call_shared_different();
					}
				}else{
					close_shared_different();
					close_select_method();
				}
				R.safari();
			};
		})(indObj);
		(function(csi_text){
			csi_text[0].style.cursor = "pointer";
			csi_text[0].onmouseover = function(){
				indObj.animate({fill:"#999"},200);
				csi_text.animate({fill:color_text_hover});
				R.safari();
			};
			csi_text[0].onmouseout = function(){
				indObj.animate({fill:color},200);
				csi_text.animate({fill:color_text});
				R.safari();
			};
			csi_text[0].onclick = function(){
				if(document.getElementById("SD_window").style.display == "none" && 
					document.getElementById("select_method").style.display == "none"){
					if(individuals[csi].fid != "0" && individuals[csi].mid != "0"){
						call_select_method();
					}else{
						call_shared_different();
					}
				}else{
					close_shared_different();
					close_select_method();
				}
				R.safari();
			};
		})(csi_text);
		initlegend();
	}
}

////////////////////////////////////////////////////////////////////
function initlegend(){
	legendset = R.set();
	legendset_share_diff = R.set();
	var lg_font = "11px \"Trebuchet MS\", Arial, sans-serif";
	var lg_x = 115;
	var lg_gap = 157;
	var lg_y = R_height-R_bottom+1;
	
	var greenobj = R.rect(lg_x,lg_y,15,15,0).attr({fill:colM,"stroke-width":0});
	var lg_green = R.text(lg_x+20,lg_y+7,"Maternal variants").attr({fill:colO,font:lg_font,"text-anchor":"start"});	
	
	blueobj = R.rect(lg_x+lg_gap*2,lg_y,15,15,0).attr({fill:colP,"stroke-width":0});
	lg_blue = R.text(lg_x+20+lg_gap*2,lg_y+7,"Paternal variants").attr({fill:colO,font:lg_font,"text-anchor":"start"});	
	
	blackobj = R.rect(lg_x+lg_gap,lg_y,15,15,0).attr({fill:colD,"stroke-width":0});
	lg_black = R.text(lg_x+20+lg_gap,lg_y+7,"De novo mutation").attr({fill:colO,font:lg_font,"text-anchor":"start"});	
	
	//var redobj = R.rect(lg_x+lg_gap*3,lg_y,15,15,0).attr({fill:colD,"stroke-width":0});
	//var lg_red = R.text(lg_x+20+lg_gap*3,lg_y+7,"De nove mutation").attr({fill:colO,font:lg_font,"text-anchor":"start"});	
	legendset.push(greenobj,lg_green,blueobj,lg_blue);
	legendset_share_diff.push(blackobj,lg_black);
	legendset.hide();
	legendset_share_diff.hide();
}

function mousePosition(ev){
	var scrollLeft = document.documentElement.scrollLeft || document.body.scrollLeft;
	var scrollTop = document.documentElement.scrollTop || document.body.scrollTop;
	return {
		x:ev.clientX + scrollLeft - document.documentElement.clientLeft,
		y:ev.clientY + scrollTop - document.documentElement.clientTop
	};
}
function mouseMove(ev){
	ev = ev || window.event;
	var mousePos = mousePosition(ev);
	restrictX = mousePos.x;
	restrictY = mousePos.y;
}
document.onmousemove = mouseMove;
document.onclick = mouseMove;


////////////////////////////////////////////////////////////////////
function show_navigator(){
	var r1 = 10;
	var r2 = 15;
	var xp = R_width + 70;
	var xm = R_width + 70;
	var y = 100;
	var delta = 4;
	var inner_l = 12;
	var inner_w = 4;

	var plus = R.set();
	plus.push(R.ellipse(xp,y-40,r2,r2).attr({fill:colO_hover,stroke:colO_hover}));
	plus.push(R.ellipse(xp,y-40,r1,r1).attr({fill:"#FFF",stroke:colO_hover}));
	plus.push(R.rect(xp-inner_w/2, y-40-inner_l/2,inner_w,inner_l,0).attr({fill:colO_hover,"stroke-width":0}));
	plus.push(R.rect(xp-inner_l/2, y-40-inner_w/2,inner_l,inner_w,0).attr({fill:colO_hover,"stroke-width":0}));
	plus.push(R.path("M"+(xp+r1-delta)+","+(y-40+r1)+
					" L"+(xp+r1-delta+r2)+","+(y-40+r1+r2)+
					" L"+(xp+r1+r2)+","+(y-40+r1-delta+r2)+
					" L"+(xp+r1)+","+(y-40+r1-delta)+"Z").attr({fill:colO_hover,"stroke-width":0}));

	(function (obj){
		for(var i = 0 ; i < obj.length ; i++){
			obj[i][0].style.cursor = "pointer";
			obj[i][0].onmouseover = function(){
				obj[1].animate({fill:colO_hover2},100);
				R.safari();
			};
			obj[i][0].onmouseout = function(){
				obj[1].animate({fill:"#FFF"},100);
				R.safari();
			};
			obj[i][0].onclick = function(){
				load_family_genome(current_chr,Math.floor(current_start+(current_end-current_start)/4),Math.ceil(current_end-(current_end-current_start)/4));
				R.safari();
			};
		}
	})(plus);

	var minus = R.set();
	minus.push(R.ellipse(xm,y,r2,r2).attr({fill:colO_hover,stroke:colO_hover}));
	minus.push(R.ellipse(xm,y,r1,r1).attr({fill:"#FFF",stroke:colO_hover}));
	minus.push(R.rect(xm-inner_l/2,y-inner_w/2,inner_l,inner_w,0).attr({fill:colO_hover,"stroke-width":0}));
	minus.push(R.path("M"+(xm+r1-delta)+","+(y+r1)+
					" L"+(xm+r1-delta+r2)+","+(y+r1+r2)+
					" L"+(xm+r1+r2)+","+(y+r1-delta+r2)+
					" L"+(xm+r1)+","+(y+r1-delta)+"Z").attr({fill:colO_hover,"stroke-width":0}));

	(function (obj){
		for(var i = 0 ; i < obj.length ; i++){
			obj[i][0].style.cursor = "pointer";
			obj[i][0].onmouseover = function(){
				obj[1].animate({fill:colO_hover2},100);
				R.safari();
			};
			obj[i][0].onmouseout = function(){
				obj[1].animate({fill:"#FFF"},100);
				R.safari();
			};
			obj[i][0].onclick = function(){
				load_family_genome(current_chr,Math.floor(current_start-(current_end-current_start)/2),Math.ceil(current_end+(current_end-current_start)/2));
				R.safari();
			};
		}
	})(minus);

	var upper = R.set();
	upper.push(R.path("M"+xm+","+(y+r2*2)+
				" L"+(xm-r1)+","+(y+r2*3)+
				" L"+(xm+r1)+","+(y+r2*3)+"Z").attr({fill:colO_hover,"stroke-width":0}));
	upper.push(R.path("M"+xm+","+(y+r2*3-delta)+
				" L"+(xm-r1)+","+(y+r2*4-delta)+
				" L"+(xm+r1)+","+(y+r2*4-delta)+"Z").attr({fill:colO_hover,"stroke-width":0}));
	(function (obj){
		for(var i = 0 ; i < obj.length ; i++){
			obj[i][0].style.cursor = "pointer";
			obj[i][0].onmouseover = function(){
				obj[0].animate({fill:colO_hover2},100);
				obj[1].animate({fill:colO_hover2},100);
				R.safari();
			};
			obj[i][0].onmouseout = function(){
				obj[0].animate({fill:colO_hover},100);
				obj[1].animate({fill:colO_hover},100);
				R.safari();
			};
			obj[i][0].onclick = function(){
				load_family_genome(current_chr,2*current_start-current_end,current_start);
				R.safari();
			};
		}
	})(upper);

	var up = R.path("M"+xm+","+(y+5+r2*4)+
				" L"+(xm-r1)+","+(y+5+r2*5)+
				" L"+(xm+r1)+","+(y+5+r2*5)+"Z").attr({fill:colO_hover,"stroke-width":0});
	(function (obj){
		obj[0].style.cursor = "pointer";
		obj[0].onmouseover = function(){
			obj.animate({fill:colO_hover2},100);
			R.safari();
		};
		obj[0].onmouseout = function(){
			obj.animate({fill:colO_hover},100);
			R.safari();
		};
		obj[0].onclick = function(){
			load_family_genome(current_chr,Math.floor(current_start-(current_end-current_start)/2),Math.ceil(current_end-(current_end-current_start)/2));
			R.safari();
		};
	})(up);

	var down = R.path("M"+xm+","+(y+15+r2*6)+
				" L"+(xm-r1)+","+(y+15+r2*5)+
				" L"+(xm+r1)+","+(y+15+r2*5)+"Z").attr({fill:colO_hover,"stroke-width":0});
	(function (obj){
		obj[0].style.cursor = "pointer";
		obj[0].onmouseover = function(){
			obj.animate({fill:colO_hover2},100);
			R.safari();
		};
		obj[0].onmouseout = function(){
			obj.animate({fill:colO_hover},100);
			R.safari();
		};
		obj[0].onclick = function(){
			load_family_genome(current_chr,Math.floor(current_start+(current_end-current_start)/2),Math.ceil(current_end+(current_end-current_start)/2));
			R.safari();
		};
	})(down);

	var downer = R.set();
	downer.push(R.path("M"+xm+","+(y+25+r2*7)+
				" L"+(xm-r1)+","+(y+25+r2*6)+
				" L"+(xm+r1)+","+(y+25+r2*6)+"Z").attr({fill:colO_hover,"stroke-width":0}));
	downer.push(R.path("M"+xm+","+(y+25+r2*8-delta)+
				" L"+(xm-r1)+","+(y+25+r2*7-delta)+
				" L"+(xm+r1)+","+(y+25+r2*7-delta)+"Z").attr({fill:colO_hover,"stroke-width":0}));
	(function (obj){
		for(var i = 0 ; i < obj.length ; i++){
			obj[i][0].style.cursor = "pointer";
			obj[i][0].onmouseover = function(){
				obj[0].animate({fill:colO_hover2},100);
				obj[1].animate({fill:colO_hover2},100);
				R.safari();
			};
			obj[i][0].onmouseout = function(){
				obj[0].animate({fill:colO_hover},100);
				obj[1].animate({fill:colO_hover},100);
				R.safari();
			};
			obj[i][0].onclick = function(){
				load_family_genome(current_chr,current_end,2*current_end-current_start);
				R.safari();
			};
		}
	})(downer);

	var lds_x = 675;
	var lds_y = 5;
	lds_btnset = R.set();
	var font_text = "15px \"Trebuchet MS\", Arial, sans-serif";
	var font_text1 = "11px \"Trebuchet MS\", Arial, sans-serif";
	//var lds_lable = R.text(190,65,"").attr({font:font_text,cursor:"pointer","font-weight":"bolder"});
	var lds_lable = R.text(lds_x,lds_y+10,"LDs").attr({font:font_text});
	//var lds_btn = R.rect(150,55,80,20,3).attr({fill:"#666",stroke:"#999",opacity:0.5,"stroke-width":0});
	
	var lds_btnrim = R.rect(lds_x+20,lds_y,51,20,2).attr({fill:"#FFF",stroke:colO_hover,opacity:1,"stroke-width":2});
	var lds_on = R.text(lds_x+32,lds_y+10,"ON").attr({font:font_text1,fill:"#FFF","font-weight":"bolder"});
	var lds_off = R.text(lds_x+59,lds_y+10,"OFF").attr({font:font_text1,fill:colO_hover,"font-weight":"bolder"});
	var lds_switch = R.rect(lds_x+23,lds_y+3,24,14,2).attr({fill:colO_hover,stroke:colO_hover,opacity:1,"stroke-width":2.5});
	var lds_btn = R.rect(lds_x+20,lds_y,50,20,2).attr({fill:"#FFF",stroke:colO_hover,opacity:0,"stroke-width":0});
	lds_btnset.push(lds_lable,lds_btnrim,lds_on,lds_off,lds_switch,lds_btn);
	lds_btnset.hide();
	if(lds_flag == 0){
		//lds_lable.attr({text:"Show LDs"});
		//lds_on.stop().animate({fill:colO_hover},200);
		lds_switch.stop().animate({transform: "t0 0"},200);
		lds_btnrim.stop().animate({fill:"#FFF"},200);
		//lds_off.stop().animate({fill: colO_hover},200);
	}else{
		//lds_lable.attr({text:"Hide LDs"});
		//lds_off.hide();
		//lds_off.stop().animate({fill:"#FFF"},200);
		lds_switch.stop().animate({transform: "t21 0"},200);
		lds_btnrim.stop().animate({fill:"#999"},200);
		//lds_on.show();
		//lds_on.stop().animate({fill:colO_hover},200);
	}
	(function (lds_btn){
		lds_btn[0].style.cursor = "pointer";
		lds_btn[0].onmouseover = function(){
		//	lds_btn.animate({fill:colO_hover2},200);
		};
		lds_btn[0].onclick = function(){
			if(lds_flag == 0){
				lds_flag = 1;
				//lds_lable.attr({text:"Hide LDs"});
				//lds_off.hide();
				lds_switch.stop().animate({transform: "t21 0"},200);
				lds_btnrim.stop().animate({fill:"#999"},200);
				//lds_on.show();
				if(csv!=-1 && variants[csv].lds != undefined){
					variants[csv].lds.show();
				}else{
					all_lds.show();
				}
			}else{
				lds_flag = 0;
				//lds_lable.attr({text:"Show LDs"});
				//lds_on.hide();
				lds_switch.stop().animate({transform: "t0 0"},200);
				lds_btnrim.stop().animate({fill:"#FFF"},200);
				//lds_off.show();
				all_lds.hide();
			}
		};
		lds_btn[0].onmouseout = function(){
		//	lds_btn.animate({fill:"#666"},200);
		};
	})(lds_btn);
	/*
	var agd_x = R_width + 70;
	var agd_y = R_height-R_bottom;
	var agd_h = 25;
	var agd_w = 60;
	//R.text(agd_x,agd_y,"Legend").attr({font:font_text,cursor:"pointer"});
	var agd_icon = R.rect(agd_x-30,agd_y-10,agd_w,agd_h,3).attr({fill:colO_hover,cursor:"pointer",stroke:colO_hover,opacity:0.5});
	(function (agd_icon){
		agd_icon[0].onmouseover = function(){
			agd_icon.animate({fill:"#999"},200);
		};
		agd_icon[0].onclick = function(){
				
		};
		agd_icon[0].onmouseout = function(){
			agd_icon.animate({fill:"#666"},200);
		};
	})(agd_icon);
	agd_icon.hide();
	*/
}

function show_denovo_muts(){
	if(req4.readyState == 4) {
		if(req4.status == 200) {
			var vsNodes = req4.responseXML.getElementsByTagName(xmlTagVariants);
			for(var i = 0 ; i < vsNodes.length ; i++){
				var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
				var vs_id = vsNodes[i].getAttribute(xmlAttributeId);
				var vPointer = 0;
				for(var j = 0 ; j < vNodes.length ; j++){
					var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
					var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
					var id = vNodes[j].getAttribute(xmlAttributeId);
					while(variants[vPointer].from < from){
						vPointer++;
					}
					if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
						if(vs_id == individuals[csi].fid){
							variants[vPointer].paternal = "N";
//							var vartable = document.getElementById("varlist_table");
//							if(vartable != undefined){
//								vartable.rows[vPointer+1].cells[5].innerHTML="N";
//							}
						}else if(vs_id == individuals[csi].mid){
							variants[vPointer].maternal = "N";
//							var vartable = document.getElementById("varlist_table");
//							if(vartable != undefined){
//								vartable.rows[vPointer+1].cells[6].innerHTML="N";
//							}
						}else if(vs_id == csi){
							change_variant_color(vPointer,colD);
						}
					}
				}
			}
		}
	}
}

function show_maternal_vars(){
	if(req5.readyState == 4) {
		if(req5.status == 200) {
			R_sremove();

			var vsNodes = req5.responseXML.getElementsByTagName(xmlTagVariants);
			for(var i = 0 ; i < vsNodes.length ; i++){
				var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
				var vs_id = vsNodes[i].getAttribute(xmlAttributeId);
				var vPointer = 0;
				for(var j = 0 ; j < vNodes.length ; j++){
					var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
					var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
					var id = vNodes[j].getAttribute(xmlAttributeId);
					while(variants[vPointer].from < from){
						vPointer++;
					}
					if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
						if(vs_id == individuals[csi].fid){
							variants[vPointer].paternal = "N";
//							var vartable = document.getElementById("varlist_table");
//							if(vartable != undefined){
//								vartable.rows[vPointer+1].cells[5].innerHTML="N";
//							}
						}else if(vs_id == individuals[csi].mid){
							variants[vPointer].maternal = "Y";
						}else if(vs_id == csi){
							change_variant_color(vPointer,colP);
						}
					}
				}
			}
		}
	}
}

function show_paternal_vars(){
	if(req6.readyState == 4) {
		if(req6.status == 200) {
			var vsNodes = req6.responseXML.getElementsByTagName(xmlTagVariants);
			for(var i = 0 ; i < vsNodes.length ; i++){
				var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
				var vs_id = vsNodes[i].getAttribute(xmlAttributeId);
				var vPointer = 0;
				for(var j = 0 ; j < vNodes.length ; j++){
					var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
					var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
					var id = vNodes[j].getAttribute(xmlAttributeId);
					while(variants[vPointer].from < from){
						vPointer++;
					}
					if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
						if(vs_id == individuals[csi].fid){
							variants[vPointer].paternal = "Y";
						}else if(vs_id == individuals[csi].mid){
							variants[vPointer].maternal = "N";
//							var vartable = document.getElementById("varlist_table");
//							if(vartable != undefined){
//								vartable.rows[vPointer+1].cells[6].innerHTML="N";
//							}
						}else if(vs_id == csi){
							change_variant_color(vPointer,colM);
						}
					}
				}
			}
		}
	}
}


function show_colorful_vars(){
	if(req5.readyState == 4) {
		if(req5.status == 200) {
			
			var bin_size = 10;
			var pixel_per_variant = 1;
			var l = R_height - R_top - R_bottom;
			var w = R_width - R_left - R_right;
			var colorbins = [];
			colorbins[0] = [];
			colorbins[1] = [];
			var int_dif_step = Math.floor((current_end - current_start)/77);
			var redbins = [];
			redbins[0] = [];
			redbins[1] = [];
			var greenbins = [];
			greenbins[0] = [];
			greenbins[1] = [];
			var bluebins = [];
			bluebins[0] = [];
			bluebins[1] = [];
			bin_set = R.set();

			recom_GB[0] = [];
			recom_GB[1] = [];
			
			var color = [];
			color[0] = colP;
			color[1] = "#83A3CD";
			color[2] = "#C1D1E6";
			color[3] = "#FFFFFF";
			color[4] = "#B3DCC4";
			color[5] = "#66BA8A";
			color[6] = colM;
			color[7] = colD;
			var colBorder = colO_hover;

			var intdif_color = [];
			intdif_color[0] = colInter;
			intdif_color[1] = "#8D6FB0";
			intdif_color[2] = "#C6B7D7";

			var help_set = R.set();
			var help = [];
			help[0] = R.path("M16,1.466C7.973,1.466,1.466,7.973,1.466,16c0,8.027,6.507,14.534,14.534,14.534c8.027,0,14.534-6.507,14.534-14.534C30.534,7.973,24.027,1.466,16,1.466z").attr({fill: "#ccc", stroke: "none",cursor     : "pointer"});
			help[1] = R.path("    M17.328,19.003v0.858h-2.707v-1.057c0-3.19,3.63-3.696,3.63-5.963c0-1.034-0.924-1.826-2.134-1.826c-1.254,0-2.354,0.924-2.354,0.924l-1.541-1.915c0,0,1.519-1.584,4.137-1.584c2.487,0,4.796,1.54     ,4.796,4.136C21.156,16.208,17.328,16.627,17.328,19.003z").attr({fill: "#000", stroke: "none",cursor: "pointer"});
			help[2] = R.path("M17.328,24.371h-2.707v-2.596h2.707V24.371z").attr({fill: "#000", stroke: "none",cursor: "pointer"});
			help_set.push(help[0],help[1],help[2]);
			help[0].stop().animate({transform: "t126 40 s0.6"});
			help[1].stop().animate({transform: "t126 40 s0.6"});
			help[2].stop().animate({transform: "t126 38 s0.8"});
			help_set.hide();
			help_set.mousedown(function(){
				window.open('tutorial.html');
			});
	
			if(!outScopeflag){
				var vsNodes = req5.responseXML.getElementsByTagName(xmlTagVariants);
				for(var i = 0 ; i < vsNodes.length ; i++){
					var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
					var vs_id = vsNodes[i].getAttribute(xmlAttributeId);
					var vPointer = 0;
					for(var j = 0 ; j < vNodes.length ; j++){
						var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
						var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
						var id = vNodes[j].getAttribute(xmlAttributeId);
						var type = vNodes[j].getAttribute(xmlAttributeType);
	
						var genotype = vNodes[j].getAttribute(xmlAttributeGenotype);
						var splitter = "|";
						if(genotype.indexOf("/")>0){
								splitter = "/";
						}
						var genotypes = genotype.split(splitter);
						
						if(id == "."){
							id = "Unnamed "+type;
						}
						while(variants[vPointer].from < from){
							vPointer++;
						}
						if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
							if(vsNodes.length == 1 && vs_id == "Intersection"){
								change_variant_color(vPointer,colInter);
								
								var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
								if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
								if(colorbins[0][Math.floor(natual_pos)] == undefined){
										colorbins[0][Math.floor(natual_pos)] = 1;
									}else{
										colorbins[0][Math.floor(natual_pos)] ++;
									}
								}
								if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
									if(colorbins[1][Math.floor(natual_pos)] == undefined){
										colorbins[1][Math.floor(natual_pos)] = 1;
									}else{
										colorbins[1][Math.floor(natual_pos)] ++;
									}
								}
	
							}else if(vsNodes.length == 1 && vs_id == "Difference"){
								change_variant_color(vPointer,colDiff);
								
								var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
								if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
									if(colorbins[0][Math.floor(natual_pos)] == undefined){
										colorbins[0][Math.floor(natual_pos)] = 1;
									}else{
										colorbins[0][Math.floor(natual_pos)] ++;
									}
								}
								if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
									if(colorbins[1][Math.floor(natual_pos)] == undefined){
										colorbins[1][Math.floor(natual_pos)] = 1;
									}else{
										colorbins[1][Math.floor(natual_pos)] ++;
									}
								}
							}else{ 
								if(vs_id == individuals[csi].fid){
									variants[vPointer].paternal = "Y";
									variants[vPointer].maternal = "N";
									change_variant_color(vPointer,colP);
									var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
									if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
										if(bluebins[0][Math.floor(natual_pos)] == undefined){
											bluebins[0][Math.floor(natual_pos)] = 1;
										}else if(genotypes[0] == "1" && genotypes[1] == "0"){
											bluebins[0][Math.floor(natual_pos)] ++;
										}else if(genotypes[0] == "0" && genotypes[1] == "1"){
											bluebins[0][Math.floor(natual_pos)] --;
										}
									}
									if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
										if(bluebins[1][Math.floor(natual_pos)] == undefined){
											bluebins[1][Math.floor(natual_pos)] = 1;
										}else if(genotypes[0] == "1" && genotypes[1] == "0"){
											bluebins[1][Math.floor(natual_pos)] ++;
										}else if(genotypes[0] == "0" && genotypes[1] == "1"){
											bluebins[1][Math.floor(natual_pos)] --;
										}
									}
								}else if(vs_id == individuals[csi].mid){
									variants[vPointer].paternal = "N"
									variants[vPointer].maternal = "Y";
									change_variant_color(vPointer,colM);
									var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
									if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
										if(greenbins[0][Math.floor(natual_pos)] == undefined){
											greenbins[0][Math.floor(natual_pos)] = 1;
										}else if(genotypes[0] == "1" && genotypes[1] == "0"){
											greenbins[0][Math.floor(natual_pos)] ++;
										}else if(genotypes[0] == "0" && genotypes[1] == "1"){
											greenbins[0][Math.floor(natual_pos)] --;
										}
									}
									if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
										if(greenbins[1][Math.floor(natual_pos)] == undefined){
											greenbins[1][Math.floor(natual_pos)] = 1;
										}else if(genotypes[0] == "1" && genotypes[1] == "0"){
											greenbins[1][Math.floor(natual_pos)] ++;
										}else if(genotypes[0] == "0" && genotypes[1] == "1"){
											greenbins[1][Math.floor(natual_pos)] --;
										}
									}
								}else if(vs_id == csi){
									variants[vPointer].paternal = "N";
									variants[vPointer].maternal = "N";
									change_variant_color(vPointer,colD);
									
									if(variants.length > l/3){
										var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*l;
										if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
											redbins[0][j] = natual_pos;
										}
										if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
											redbins[1][j] = natual_pos;
										}
									}
									/*var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
									if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
										if(redbins[0][Math.floor(natual_pos)] == undefined){
											redbins[0][Math.floor(natual_pos)] = 1;
										}else{
											redbins[0][Math.floor(natual_pos)] ++;
										}
									}
									if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
										if(redbins[1][Math.floor(natual_pos)] == undefined){
											redbins[1][Math.floor(natual_pos)] = 1;
										}else{
											redbins[1][Math.floor(natual_pos)] ++;
										}
									}*/
								}
							}
						}
					}
				}
			}else{
				var vlNodes = req5.responseXML.getElementsByTagName(xmlTagValues);
				for(var i=0; i<vlNodes.length; i++){
					if(vlNodes[i].getAttribute(xmlAttributeId) == "IntersectionL" || vlNodes[i].getAttribute(xmlAttributeId) == "DifferenceL"){
						colorbins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
						int_dif_step = vlNodes[i].getElementsByTagName(xmlTagStep)[0].firstChild.nodeValue;
					}else if(vlNodes[i].getAttribute(xmlAttributeId) == "IntersectionR" || vlNodes[i].getAttribute(xmlAttributeId) == "DifferenceR"){
						colorbins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
						int_dif_step = vlNodes[i].getElementsByTagName(xmlTagStep)[0].firstChild.nodeValue;
					}else{
						switch(vlNodes[i].getAttribute(xmlAttributeId)){
							case (families[csi].fid + "L"):
								bluebins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
							case (families[csi].fid + "R"):
								bluebins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
							case (families[csi].mid + "L"):
								greenbins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
							case (families[csi].mid + "R"):
								greenbins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
						}
					}
				}

				/*var vsNodes = req5.responseXML.getElementsByTagName(xmlTagVariants);
				for(var i=0; i<vsNodes.length; i++){
					var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
					for(var j = 0 ; j < vNodes.length ; j++){
						var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
						var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;

						var genotype = vNodes[j].getAttribute(xmlAttributeGenotype);
						var splitter = "|";
						if(genotype.indexOf("/")>0){
							splitter = "/";
						}
						var genotypes = genotype.split(splitter);

						if(from >= current_start && to <= current_end){
							var natual_pos = ((parseInt(to)+parseInt(from))/2 - current_start)/(current_end - current_start + 1)*l;
							if(genotypes[0] != undefined && genotypes[0] != "0"){
								redbins[0][j] = natual_pos;
								//bin_set.push(R.rect(R_left+w/2-40+6,natual_pos+R_top,bin_size,2,0).attr({fill:colD,stroke:colD,"stroke-width":1}).toFront());
							}
							if(genotypes[1] != undefined && genotypes[1] != "0"){
								redbins[1][j] = natual_pos;
								//bin_set.push(R.rect(R_left+w/2+40-16,natual_pos+R_top,bin_size,2,0).attr({fill:colD,stroke:colD,"stroke-width":1}).toFront());
							}
						}
					}
				}*/
			}
			if(outScopeflag || variants.length >= l/3){
				var ratio = 0;
				var re_i = 0;

				if(greenbins[0].length != 0 || greenbins[1].length != 0 || bluebins[0].length != 0 || bluebins[1].length != 0){
					help_set.show();
					bin_set.push(R.rect(R_left+w/2-40+6-100,R_top-30,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
					bin_set.push(R.rect(R_left+w/2-40+6-90,R_top-30,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
					bin_set.push(R.rect(R_left+w/2-40+6-80,R_top-30,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
					bin_set.push(R.rect(R_left+w/2-40+6-70,R_top-30,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
				}
				if(colorbins[0].length != 0){
					help_set.show();
					bin_set.push(R.rect(R_left+w/2-40+6-100,R_top-30,10,bin_size,0).attr({fill:intdif_color[0],stroke:colBorder,"stroke-width":1}));
					bin_set.push(R.rect(R_left+w/2-40+6-90,R_top-30,10,bin_size,0).attr({fill:intdif_color[1],stroke:colBorder,"stroke-width":1}));
					bin_set.push(R.rect(R_left+w/2-40+6-80,R_top-30,10,bin_size,0).attr({fill:intdif_color[2],stroke:colBorder,"stroke-width":1}));
				}

				var if_firstB = [];
				var if_firstG = [];
				for(var i=0; i<bins[0].length; i++){
					if(colorbins[0][i] > 0){
						ratio = colorbins[0][i] / parseInt(int_dif_step);
						if(10000*ratio > 1){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[0],stroke:colBorder,"stroke-width":1}));
						}else if(100000*ratio > 1){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[1],stroke:colBorder,"stroke-width":1}));
						}else{
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[2],stroke:colBorder,"stroke-width":1}));
						}
					}
					if(greenbins[0][i] == undefined){
						greenbins[0][i] = 0;
					}
					if(bluebins[0][i] == undefined){
						bluebins[0][i] = 0;
					}
					if(Math.abs(bluebins[0][i]) > 0 || Math.abs(greenbins[0][i]) > 0){
						if(Math.abs(bluebins[0][i]) >= Math.abs(greenbins[0][i])){
							if(if_firstB[0] == undefined){
								if(bluebins[0][i] > 0){
									if_firstB[0] = 1;
								}else{
									if_firstB[0] = -1;
								}
							}	
							if(bluebins[0][i] > 0){
								if(if_firstB[0] > 0){
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 1;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 2;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstB[0] > 0){
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 2;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 1;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}
							}
						}else if(Math.abs(bluebins[0][i]) < Math.abs(greenbins[0][i])){
							if(if_firstG[0] == undefined){
								if(greenbins[0][i] > 0){
									if_firstG[0] = 1;
								}else{
									if_firstG[0] = -1;
								}
							}
							if(greenbins[0][i] > 0){
								if(if_firstG[0] > 0){
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 4;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 3;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstG[0] > 0){
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 3;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[0][re_i]={};
									recom_GB[0][re_i].value = 4;
									recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}
							}
						}
					}

					/*
					if(bluebins[0][i] > 0 || greenbins[0][i] > 0){
						ratio = greenbins[0][i] / (parseInt(greenbins[0][i],10)+parseInt(bluebins[0][i],10));
						if(7*ratio < 1){
							recom_GB[0][re_i]={};
							recom_GB[0][re_i].value = 1;
							recom_GB[0][re_i].pos = i;
							re_i = re_i+1;
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 2){
							recom_GB[0][re_i]={};
							recom_GB[0][re_i].value = 1;
							recom_GB[0][re_i].pos = i;
							re_i = re_i+1;
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 3){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[2],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 4){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[3],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 5){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[4],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 6){
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
							recom_GB[0][re_i]={};
							recom_GB[0][re_i].value = 2;
							recom_GB[0][re_i].pos = i;
							re_i = re_i+1;
						}else{
							bin_set.push(R.rect(R_left+w/2-40+6,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
							recom_GB[0][re_i]={};
							recom_GB[0][re_i].value = 2;
							recom_GB[0][re_i].pos = i;
							re_i = re_i+1;
						}
					}
					*/
					//if(redbins[0][i] != undefined){
					//	bin_set.push(R.rect(R_left+w/2-40+18,i*bin_size+R_top,2,bin_size,0).attr({fill:color[7],stroke:"#FFFFFF","stroke-width":1}));
					//}
					
				}

				re_i = 0;
				for(var i=0; i<bins[1].length; i++){
					if(colorbins[1][i] > 0){
						ratio = colorbins[1][i] / parseInt(int_dif_step);
						if(10000*ratio > 1){
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[0],stroke:colBorder,"stroke-width":1}));
						}else if(100000*ratio > 1){
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[1],stroke:colBorder,"stroke-width":1}));
						}else{
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:intdif_color[2],stroke:colBorder,"stroke-width":1}));
						}
					}
					if(greenbins[1][i] == undefined){
						greenbins[1][i] = 0;
					}
					if(bluebins[1][i] == undefined){
						bluebins[1][i] = 0;
					}
					if(Math.abs(bluebins[1][i]) > 0 || Math.abs(greenbins[1][i]) > 0){
						if(Math.abs(bluebins[1][i]) >= Math.abs(greenbins[1][i])){
							if(if_firstB[0] == undefined){
								if(bluebins[1][i] > 0){
									if_firstB[0] = 1;
								}else{
									if_firstB[0] = -1;
								}
							}	
							if(bluebins[1][i] > 0){
								if(if_firstB[0] > 0){
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 1;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 2;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstB[0] > 0){
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 2;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 1;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}
							}
						}else if(Math.abs(bluebins[1][i]) < Math.abs(greenbins[1][i])){
							if(if_firstG[0] == undefined){
								if(greenbins[1][i] > 0){
									if_firstG[0] = 1;
								}else{
									if_firstG[0] = -1;
								}
							}
							if(greenbins[1][i] > 0){
								if(if_firstG[0] > 0){
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 4;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 3;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstG[0] > 0){
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 3;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}else{
									recom_GB[1][re_i]={};
									recom_GB[1][re_i].value = 4;
									recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}
							}
						}
					}
					/*
					if(bluebins[1][i] > 0 || greenbins[1][i] > 0){
						ratio = greenbins[1][i] / (parseInt(greenbins[1][i],10)+parseInt(bluebins[1][i],10));
						if(7*ratio < 1){
							recom_GB[1][re_i]={};
							recom_GB[1][re_i].value = 1;
							recom_GB[1][re_i].pos = i;
							re_i = re_i+1;
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 2){
							recom_GB[1][re_i]={};
							recom_GB[1][re_i].value = 1;
							recom_GB[1][re_i].pos = i;
							re_i = re_i+1;
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 3){
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[2],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 4){
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[3],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 5){
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[4],stroke:colBorder,"stroke-width":1}));
						}else if(7*ratio < 6){
							recom_GB[1][re_i]={};
							recom_GB[1][re_i].value = 2;
							recom_GB[1][re_i].pos = i;
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
							re_i = re_i+1;
						}else{
							recom_GB[1][re_i]={};
							recom_GB[1][re_i].value = 2;
							recom_GB[1][re_i].pos = i;
							bin_set.push(R.rect(R_left+w/2+40-16,i*bin_size+R_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
							re_i = re_i+1;
						}
					}
					*/
					//if(redbins[1][i] != undefined){
					//	bin_set.push(R.rect(R_left+w/2+40-20,i*bin_size+R_top,2,bin_size,0).attr({fill:color[7],stroke:"#FFFFFF","stroke-width":1}));
					//}
				}
				for(var j=0; j<redbins[0].length; j++){
					if(redbins[0][j] != undefined){
						bin_set.push(R.rect(R_left+w/2-40+6,redbins[0][j]+R_top,bin_size,2,0).attr({fill:colD,stroke:colD,"stroke-width":1}).toFront());
					}
				}
				for(var j=0; j<redbins[1].length; j++){
					if(redbins[1][j] != undefined){
						bin_set.push(R.rect(R_left+w/2+40-16,redbins[1][j]+R_top,bin_size,2,0).attr({fill:colD,stroke:colD,"stroke-width":1}).toFront());
					}
				}

				var chr_left;
				var chr_right;
				if(current_end - current_start <= 20000000){
					chr_left = newif_recombination(recom_GB[0],5);
					chr_right = newif_recombination(recom_GB[1],5);
				}else{
					chr_left = newif_recombination(recom_GB[0],3);
					chr_right = newif_recombination(recom_GB[1],3);
				}
				var rcbicon1 = R.path("M26.711,14.086L16.914,4.29c-0.778-0.778-2.051-0.778-2.829,0L4.29,14.086c-0.778,0.778-0.778,2.05,0,2.829l9.796,9.796c0.778,0.777,2.051,0.777,2.829,0l9.797-9.797C27.488,16.136,27.488,14.864,26.711,14.086z").attr({fill: "#ccc", stroke: "none",cursor: "pointer"});
				var rcbicon2 = R.path("M14.702,8.981c0.22-0.238,0.501-0.357,0.844-0.357s0.624,0.118,0.844,0.353c0.221,0.235,0.33,0.531,0.33,0.885c0,0.306-0.101,1.333-0.303,3.082c-0.201,1.749-0.379,3.439-0.531,5.072H15.17c-0.135-1.633-0.301-3.323-0.5-5.072c-0.198-1.749-0.298-2.776-0.298-3.082C14.372,9.513,14.482,9.22,14.702,8.981z").attr({fill: "#000", stroke: "none",cursor: "pointer"});
				var rcbicon3 = R.path("M16.431,21.799c-0.247,0.241-0.542,0.362-0.885,0.362s-0.638-0.121-0.885-0.362c-0.248-0.241-0.372-0.533-0.372-0.876s0.124-0.638,0.372-0.885c0.247-0.248,0.542-0.372,0.885-0.372s0.638,0.124,0.885,0.372c0.248,0.247,0.372,0.542,0.372,0.885S16.679,21.558,16.431,21.799z").attr({fill: "#000", stroke: "none",cursor: "pointer"});
				rcbicon1.hide();
				rcbicon2.hide();
				rcbicon3.hide();
				recom_set = R.set();
				var icon_set = R.set();
				recom_set.push(rcbicon1,rcbicon2,rcbicon3);
				if(chr_left.length != 0){
					for(var i=0; i<chr_left.length; i++){
						var clone1 = rcbicon1.clone();
						var clone2 = rcbicon2.clone();
						var clone3 = rcbicon3.clone();
						clone1[0].onclick = function(){
							call_siblings_window();
						};
						clone2[0].onclick = function(){
							call_siblings_window();
						};
						clone3[0].onclick = function(){
							call_siblings_window();
						};
						clone1.animate({transform:"t"+(R_left+w/2-40+16)+" "+(R_top+chr_left[i].pos*bin_size-6)});
						clone2.animate({transform:"t"+(R_left+w/2-40+16)+" "+(R_top+chr_left[i].pos*bin_size-6)});
						clone3.animate({transform:"t"+(R_left+w/2-40+16)+" "+(R_top+chr_left[i].pos*bin_size-6)});
						recom_set.push(clone1,clone2,clone3);
						icon_set.push(clone2,clone3);
					}
				}
				if(chr_right.length != 0){
					for(var i=0; i<chr_right.length; i++){
						var clone1 = rcbicon1.clone();
						var clone2 = rcbicon2.clone();
						var clone3 = rcbicon3.clone();
						clone1[0].onclick = function(){
							call_siblings_window();
						};
						clone2[0].onclick = function(){
							call_siblings_window();
						};
						clone3[0].onclick = function(){
							call_siblings_window();
						};
						clone1.animate({transform:"t"+(R_left+w/2+40-46)+" "+(R_top+chr_right[i].pos*bin_size-6)});
						clone2.animate({transform:"t"+(R_left+w/2+40-46)+" "+(R_top+chr_right[i].pos*bin_size-6)});
						clone3.animate({transform:"t"+(R_left+w/2+40-46)+" "+(R_top+chr_right[i].pos*bin_size-6)});
						recom_set.push(clone1,clone2,clone3);
						icon_set.push(clone2,clone3);
					}
				}

				var anim_0 = Raphael.animation({
					100:{
						fill:"#000"
					},
					50:{
						fill:"#ccc"
					}},1000);
				icon_set.animate(anim_0.repeat("Infinity"));
				

				/*
				if(chr_left.length != 0 && chr_right.length != 0 && chr_left.length == chr_right.length && chr_left[0].value != chr_right[0].value){
					var leftcolor;
					var rightcolor;
					var leftcolor_hover;
					var rightcolor_hover;
					var leftchar = "";
					var rightchar = "";
					var l = R_height - R_top - R_bottom;
					var w = R_width - R_left - R_right;
					if(chr_left[0].value == 1){
						leftcolor = colP;
						rightcolor = colM;
						leftcolor_hover = colP_hover;
						rightcolor_hover = colM_hover;
						leftchar = "M"+ (R_left+w/2+30) +","+ 15 + " L"+ (R_left+w/2-10) +","+ R_top;
						rightchar = "M"+(R_left+w/2-30)+","+ 15 + " L"+ (R_left+w/2+10) +","+ R_top;
					}else{
						leftcolor = colM;
						rightcolor = colP;
						leftcolor_hover = colM_hover;
						rightcolor_hover = colP_hover;
						leftchar = "M"+(R_left+w/2-30)+","+ 15 +" L"+ (R_left+w/2-10) +","+ R_top;
						rightchar = "M"+ (R_left+w/2+30) +","+ 15 + " L"+ (R_left+w/2+10) +","+ R_top;
					}
					for(var i=0; i<chr_left.length; i++){
						var gap_left = chr_left[i].nextpos - chr_left[i].pos;
						var gap_right = chr_right[i].nextpos - chr_right[i].pos;
						var recomb_pos = -1;
						var recomb_nextpos = -1;
						if(gap_left > gap_right){
							recomb_pos = chr_right[i].pos;
							recomb_nextpos = chr_right[i].nextpos;
						}else{
							recomb_pos = chr_left[i].pos;
							recomb_nextpos = chr_left[i].nextpos;
						}
						if(i%2==0){
							leftchar = leftchar + " L"+ (R_left+w/2-10) +","+ (R_top+recomb_pos*bin_size) + 
										" L"+ (R_left+w/2+10) +","+ (R_top+recomb_nextpos*bin_size);
							rightchar = rightchar + " L"+ (R_left+w/2+10) +","+ (R_top+recomb_pos*bin_size) +
										" L"+ (R_left+w/2-10) +","+ (R_top+recomb_nextpos*bin_size);
						}else{
							leftchar = leftchar + " L"+ (R_left+w/2+10) +","+ (R_top+recomb_pos*bin_size) +
										" L"+ (R_left+w/2-10) +","+ (R_top+recomb_nextpos*bin_size);
							rightchar = rightchar + " L"+ (R_left+w/2-10) +","+ (R_top+recomb_pos*bin_size) +
										" L"+ (R_left+w/2+10) +","+ (R_top+recomb_nextpos*bin_size);
						}
					}
					if(chr_left[0].value == 1){
						if(chr_left[chr_left.length-1].value == 2){
							leftchar = leftchar + " L"+ (R_left+w/2-10) +","+ (R_top+l);
							rightchar = rightchar + " L"+ (R_left+w/2+10) +","+ (R_top+l);
						}else{
							leftchar = leftchar + " L"+ (R_left+w/2+10) +","+ (R_top+l);
							rightchar = rightchar + " L"+ (R_left+w/2-10) +","+ (R_top+l);
						}
					}else{
						if(chr_left[chr_left.length-1].value == 1){
							leftchar = leftchar + " L"+ (R_left+w/2-10) +","+ (R_top+l);
							rightchar = rightchar + " L"+ (R_left+w/2+10) +","+ (R_top+l);
						}else{
							leftchar = leftchar + " L"+ (R_left+w/2+10) +","+ (R_top+l);
							rightchar = rightchar + " L"+ (R_left+w/2-10) +","+ (R_top+l);
						}
					}
					recom_set = R.set();
					var colorline = [];
					colorline[0] = R.path(leftchar).attr({stroke:leftcolor,"stroke-width":6,opacity:1,cursor:"pointer"});
					colorline[1] = R.path(rightchar).attr({stroke:rightcolor,"stroke-width":6,opacity:1,cursor:"pointer"});
					
					
					
					
					/*
					var chr_left = if_recombination(recom_GB[0]);
					var chr_right = if_recombination(recom_GB[1]);

					var recomb_pos = -1;
					var recomb_uppos = -1;
					var gap_left = recomb_loaction(recom_GB[0]);
					var gap_right = recomb_loaction(recom_GB[1]); 
					if(gap_left.gap > gap_right.gap){
						recomb_pos = gap_right.pos;
						recomb_uppos = gap_right.uppos;
					}else{
						recomb_pos = gap_left.pos;
						recomb_uppos = gap_left.uppos;
					}
					*/
					
					/*
					recom_set = R.set();
					var l = R_height - R_top - R_bottom;
					var w = R_width - R_left - R_right;
					var colorline = [];
					if(chr_left==1 && chr_right==2){
						colorline[0] = R.path("M"+ (R_left+w/2+30) +","+ 15 +" L"+ (R_left+w/2-10) +","+ R_top +
								" L"+ (R_left+w/2-10) +","+ (R_top+recomb_uppos*bin_size) +" L"+ (R_left+w/2+10) +","+ (R_top+bin_size+recomb_pos*bin_size) +
								" L"+ (R_left+w/2+10) +","+ (R_top+l)
								).attr({stroke:colP,"stroke-width":6,opacity:1,cursor:"pointer"});
						colorline[1] = R.path("M"+(R_left+w/2-30)+","+ 15 +" L"+ (R_left+w/2+10) +","+ R_top +
								" L"+ (R_left+w/2+10) +","+ (R_top+recomb_uppos*bin_size) +" L"+ (R_left+w/2-10) +","+ (R_top+bin_size+recomb_pos*bin_size) +
								" L"+ (R_left+w/2-10) +","+ (R_top+l)
								).attr({stroke:colM,"stroke-width":6,opacity:1,cursor:"pointer"});
					}else if(chr_left==2 && chr_right==1){
						colorline[0] = R.path("M"+ (R_left+w/2-30) +","+ 15 +" L"+ (R_left+w/2-10) +","+ R_top +
								" L"+ (R_left+w/2-10) +","+ (R_top+recomb_uppos*bin_size) +" L"+ (R_left+w/2+10) +","+ (R_top+bin_size+recomb_pos*bin_size) +
								" L"+ (R_left+w/2+10) +","+ (R_top+l)
								).attr({stroke:colP,"stroke-width":6,opacity:1,cursor:"pointer"});
						colorline[1] = R.path("M"+(R_left+w/2+30)+","+ 15 +" L"+ (R_left+w/2+10) +","+ R_top +
								" L"+ (R_left+w/2+10) +","+ (R_top+recomb_uppos*bin_size) +" L"+ (R_left+w/2-10) +","+ (R_top+bin_size+recomb_pos*bin_size) +
								" L"+ (R_left+w/2-10) +","+ (R_top+l)
								).attr({stroke:colM,"stroke-width":6,opacity:1,cursor:"pointer"});
					}
					*/
				/*	
				 18yy	var anim_0 = Raphael.animation({
						50:{
							"stroke-width":8,
							stroke:leftcolor_hover
						},
						100:{
							"stroke-width":6,
							stroke:leftcolor
						}},2000);
					var anim_1 = Raphael.animation({
						50:{
							"stroke-width":8,
							stroke:rightcolor_hover
						},
						100:{
							"stroke-width":6,
							stroke:rightcolor
						}},2000);
					colorline[0].animate(anim_0.repeat("Infinity"));
					colorline[1].animate(anim_1.repeat("Infinity"));
					recom_set.push(colorline[0],colorline[1]);
					recom_set.mouseover(function(){
						colorline[0].animate({stroke:leftcolor_hover},200);
						colorline[1].animate({stroke:rightcolor_hover},200);
					});
					recom_set.mouseout(function(){
						colorline[0].animate({stroke:leftcolor},200);
						colorline[1].animate({stroke:rightcolor},200);
					});
					colorline[0][0].onclick = function(){
						call_siblings_window();
					};
					colorline[1][0].onclick = function(){
						call_siblings_window();
					};
				}
				*/
			}else{
				help_set.hide();
			}
			
			if(compare_method == "trioAnalysis"){
				legendset.show();
				blueobj.attr({fill:colP});
				lg_blue.attr({text:"Paternal variants"});
				lg_black.attr({text:"De nove mutation"});
				blackobj.attr({fill:colD});
				legendset_share_diff.animate({transform:"t0 0"})
				legendset_share_diff.show();
			}
			R_sremove();
			warning_set.hide();
		}
	}
}
function clear_colorful_vars(){
	for(var i = 0 ; i < variants.length ; i++){
		change_variant_color(i,colO);
	}
	if(bin_set != null){
		bin_set.remove();
	}
	if(recom_set != null){
		recom_set.remove();
	}
}

function family_switch_compareindiv(){
	parents_area.hide();
}

function highlight_vars(){
	var font_size = "14px \"Trebuchet MS\", Arial, sans-serif";
	var w = R_width - R_left - R_right;
	var len = 0;
	for(id in compared_individuals){
		if(id == csi){
			delete compared_individuals[id];
		}else{
			len++;
		}
	}
	clear_colorful_vars();
	close_shared_different();
	//family_switch_compareindiv();
	parents_area.hide();
	compared_indiv_area.remove();
	if((compare_method == "getDifference" || compare_method == "getIntersection") 
		&& len > 0 && len <= 9 && individuals[csi] != undefined){
		var sets5 = "";
		for(var id in compared_individuals){
			sets5 += id + ":";
		}

		var temp_str = "";
		if(len < 4){
			for(var i=0; i<len; i++){
				if(i == len-1){
					temp_str += sets5.split(":")[i];
				}else{
					temp_str += sets5.split(":")[i] + ",";
				}
			}
		}else{
			temp_str = sets5.split(":")[0] + "," + sets5.split(":")[1] + "," + sets5.split(":")[2] + "..."
		}
		compared_indiv_area.push(R.text(R_left+w/2+2,R_top/2-5,temp_str).attr({font:font_size,"text-anchor":"middle"}));


		if(compare_method == "getDifference"){
			sets5 = sets5.substring(0,sets5.length-1) + ",";
			//csi_text.attr({fill:colDiff});
			compared_indiv_area.push(R.text(R_left+w/2+2,R_top/2-25,"Different").attr({font:font_size,"text-anchor":"middle",fill:colDiff}));
			legendset.hide();
			blackobj.attr({fill:colInter});
			lg_black.attr({text:"Unique with all selected individuals"});
			//lg_black.attr({text:"Shared variants"});
			//blackobj.attr({fill:colO});
			legendset_share_diff.animate({transform:"t-90 0"})
			legendset_share_diff.show();
		}else{
			//compared_indiv_area.attr({fill:colInter});
			//csi_text.attr({fill:colInter});
			compared_indiv_area.push(R.text(R_left+w/2+2,R_top/2-25,"Share").attr({font:font_size,"text-anchor":"middle",fill:colInter}));
			legendset.hide();
			blackobj.attr({fill:colInter});
			lg_black.attr({text:"Shared with all selected individuals"});
			//lg_black.attr({text:"different variants"});
			//blackobj.attr({fill:colO});
			legendset_share_diff.animate({transform:"t-90 0"})
			legendset_share_diff.show();
		}
		sets5 += csi;
		var form5 = new FormData();
		form5.append("sets",sets5);
		form5.append("enctype","multipart/form-data");
		req5.onreadystatechange = show_colorful_vars;
		querry = "action="+compare_method+"&tracks="+trackname+"&chr="+current_chr+"&start="+current_start+"&end="+current_end;
		req5.open("POST","servlet/test.do?"+querry,true);
		req5.send(form5);
	}else{
		R_sremove();
	}
}
function trioAnalysis(){
	if(compare_method == "trioAnalysis" 
		&& individuals[csi].fid != "0" && individuals[csi].mid != "0"){
		clear_colorful_vars();
		req5.onreadystatechange = show_colorful_vars;
		querry = "action=trioAnalysis&tracks="+trackname+"&chr="+current_chr+"&start="+current_start+"&end="+current_end+"&id="+csi;
		req5.open("GET","servlet/test.do?"+querry,true);
		req5.send();
	}else{
		R_sremove();
		warning_set.hide();
	}
}
function call_compared_parents(){
//	parents_area.show();
	compare_method = "trioAnalysis";
	parents_area.show();
	compared_indiv_area.remove();
	trioAnalysis();
	close_shared_different();
	close_select_method();
}

function show_vars(){
	if(req3.readyState == 4) {
		if(req3.status == 200) {
			gdf_list = [];	
			if(!outScopeflag){
				init_individual_vars();
			}
			else{
				init_genes_omim();
			}
			list_variants();
			
			plot_variants();
			plot_genes();
			
			/*
			var sets4 = individuals[csi].fid
					+":"+individuals[csi].mid
					+","+csi;
			var sets5 = individuals[csi].fid
					+","+individuals[csi].mid
					+":"+csi;
			var sets6 = individuals[csi].mid
					+","+individuals[csi].fid
					+":"+csi;
			*/
			if(!outScopeflag && variants.length < (R_height - R_top - R_bottom)/15){
				req7.onreadystatechange = show_ld_relationships;
				querry = "action=getLD&tracks="+trackname+"&start="+current_start+"&end="+current_end+"&id="+csi;
				req7.open("GET","servlet/test.do?"+querry,true);
				req7.send();
			}
			R_sremove();
			if(compare_method == "trioAnalysis"){
				R_sremove = spinner(R);
				trioAnalysis();
			}else if(compare_method == "getIntersection" || compare_method == "getDifference"){
				R_sremove = spinner(R);
				highlight_vars();
			}
			/*
			var form5 = new FormData();
			form5.append("sets",sets5);
			form5.append("enctype","multipart/form-data");
			req5.onreadystatechange = show_maternal_vars;
			querry = "action=getDifference&tracks="+trackname+"&chr="+current_chr+"&start="+current_start+"&end="+current_end;
			req5.open("POST","servlet/test.do?"+querry,true);
			req5.send(form5);
		
			var form6 = new FormData();
			form6.append("sets",sets6);
			req6.onreadystatechange = show_paternal_vars;
			querry = "action=getDifference&tracks="+trackname+"&chr="+current_chr+"&start="+current_start+"&end="+current_end;
			req6.open("POST","servlet/test.do?"+querry,true);
			req6.send(form6);

			var form4 = new FormData();
			form4.append("sets",sets4);
			form4.append("enctype","multipart/form-data");
			req4.onreadystatechange = show_denovo_muts;
			querry = "action=getDifference&tracks="+trackname+"&chr="+current_chr+"&start="+current_start+"&end="+current_end;
			req4.open("POST","servlet/test.do?"+querry,true);
			req4.send(form4);
			*/
		}else{
			var pattern = /<.*?>/g;
			XMLHttpReq12.open("GET","servlet/test.do?action=getCheck&tracks=_"+initPvar_superid,false);
			XMLHttpReq12.send(null);
			var err_text=XMLHttpReq12.responseText.replace(pattern,"");
			if(err_text==null||err_text==""){
				alert("An error occurred while executing the query!\nPlease check the data format and data index.");
			}else{
				alert(err_text);
			}
		}
	}
}
function show_ld_relationships(){
	if(req7.readyState == 4) {
		if(req7.status == 200) {
			var pattern = /<.*?>/g;
			var lds = req7.responseText.replace(pattern,"").split(",");
			var l = R_height - R_top - R_bottom;
			var w = R_width - R_left - R_right;
			for(var i = 0 ; i < lds.length ; i++){
				var ld = lds[i].split(";");
				if(variants_byid[ld[0]] != undefined && variants_byid[ld[1]] != undefined){
					var var1_pos = ((variants[variants_byid[ld[0]]].to+variants[variants_byid[ld[0]]].from)/2 - current_start + 0.5)/(current_end - current_start + 1)*l+R_top;
					var var2_pos = ((variants[variants_byid[ld[1]]].to+variants[variants_byid[ld[1]]].from)/2 - current_start + 0.5)/(current_end - current_start + 1)*l+R_top;
					var x = R_left+w/2-30;
					var direction = 0;
					if(ld[2] == 1){
						x = R_left+w/2+30;
						direction = 1;
					}
					var obj_temp = R.path("M"+x+","+var1_pos+" A"+(var2_pos-var1_pos)+","+(var2_pos-var1_pos)+" 0 0,"+direction+" "+x+","+var2_pos).attr({stroke:Raphael.getRGB("rgb(255,"+(1-ld[3])*255+","+(1-ld[3])*255+")").hex});
					if(variants[variants_byid[ld[0]]].lds == undefined){
						variants[variants_byid[ld[0]]].lds = R.set();
					}
					if(variants[variants_byid[ld[1]]].lds == undefined){
						variants[variants_byid[ld[1]]].lds = R.set();
					}
					variants[variants_byid[ld[0]]].lds.push(obj_temp);
					variants[variants_byid[ld[1]]].lds.push(obj_temp);
//					obj_temp.hide(); // for hide all lds in start, and show lds when select an individual
					all_lds.push(obj_temp); //for show all lds in start, and hide inconnected lds when select an individual
					if(lds_flag == 0){
						all_lds.hide();
					}else{
						all_lds.show();
					}
				}
			}
		}
	}
}
function list_variants(){
	document.getElementById("varlist").innerHTML= "";
	if(!outScopeflag){
		var temp = document.createElement("table");
		document.getElementById("varlist").appendChild(temp);
		temp.className = "listt_table";
		temp.id = "varlist_table";
		var temp_tr = temp.insertRow(-1);
		temp_tr.innerHTML = 
			"<th>ID</th>"+
			"<th>Type</th>"+
			"<th>Pos</th>"+
			"<th>Alt</th>"+
			"<th>Genotype</th>"+
	//		"<th>Paternal</th>"+
	//		"<th>Maternal</th>"+
			"<th>AA change</th>"+
			"<th>Selected</th>"+
			"<th>DBSNP</th>";
		for(var i=0 ; i<variants.length ; i++){ 
			temp_tr = temp.insertRow(-1);
			temp_tr.style.background = "#CCC";
			(function(i){
				var temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.innerHTML = variants[i].id;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.innerHTML = variants[i].type;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.innerHTML = variants[i].chr+":"+variants[i].from+"-"+variants[i].to;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.style.width = "100px";
				temp_td.style.wordBreak = "break-all";
				temp_td.innerHTML = variants[i].letter;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.innerHTML = variants[i].genotype;
		//		temp_td = temp_tr.insertCell(-1);
		//		temp_td.innerHTML = variants[i].paternal;
		//		temp_td = temp_tr.insertCell(-1);
		//		temp_td.innerHTML = variants[i].maternal;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				temp_td.innerHTML = variants[i].functional;
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";

				var radioObj = document.createElement("input");
				temp_td.appendChild(radioObj);
				radioObj.type = "radio";
				radioObj.name = "varlist_select";
				radioObj.id = i + "__varlist_select";
				if(variants[i].id == csv){
					radioObj.checked = true;
				}
				radioObj.onclick = function(event){
					var target = event.target || event.srcElement;
					var idx = target.getAttribute("id").split("__")[0];
					select_a_variant(idx);
				};
				temp_td = temp_tr.insertCell(-1);
				temp_td.align = "center";
				if(variants[i].dd != null){
					temp_td.innerHTML = variants[i].dd;
					temp_td.style.textDecoration="underline";
					temp_td.style.cursor="pointer";
					temp_td.onmouseover = function(){
						temp_td.style.color="#FF0";
					};
					temp_td.onmouseout = function(){
						temp_td.style.color="#000";
					};
					temp_td.onclick = function(){
						window.open(refSNPurl+variants[i].dd);
					};
				}else{
					temp_td.innerHTML = "--";
				}
			})(i);
		}
	}else{
		document.getElementById("varlist").innerHTML= "The number of variants is too huge to show! The list of variants is available within the range of 300000"
	}
}

function plot_variants(){
	var l = R_height - R_top - R_bottom;
	var w = R_width - R_left - R_right;
	var font_height = 15;
	var text_length = 100;
	var bin_size = 10;
	var max_bar_length = 100;
	var font_size2_text = "12px \"Trebuchet MS\", Arial, sans-serif";
	var pixel_per_variant = 1;

	if(outScopeflag){
		var vlNodes = req3.responseXML.getElementsByTagName(xmlTagValues);
		var num_temp1 = 0;
		var num_temp2 = 0;
		for(var i=0; i<vlNodes.length; i++){
			if(vlNodes[i].getAttribute(xmlAttributeId).charAt(vlNodes[i].getAttribute(xmlAttributeId).length-1)=="L"){
					//alert("L");
				var nums_temp = [];
				nums_temp = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
				for(var j = 0 ; j < nums_temp.length ; j++){
					num_temp1 = num_temp1 + parseInt(nums_temp[j]);
				}
			}else{
				var nums_temp = [];
				nums_temp = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
				for(var j = 0 ; j < nums_temp.length ; j++){
					num_temp2 = num_temp2 + parseInt(nums_temp[j]);
				}
			}
		}
		if(num_temp1 > 2* num_temp2 || num_temp2 > 2* num_temp1){
			R.rect(R_left+w/2-5,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
		}
		else{
			R.rect(R_left+w/2-15,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
			R.rect(R_left+w/2+5,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
		}
	}
	else if(variants[0] != undefined && variants[0].genotype.indexOf("|")>=0){
		R.rect(R_left+w/2-15,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
		R.rect(R_left+w/2+5,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
	} else {
		R.rect(R_left+w/2-5,R_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
	}
	if(outScopeflag || variants.length >= l/3){
		lds_btnset.hide();
		bins[0] = [];
		bins[1] = [];
		var maxbins = [];
		maxbins[0] = 0;
		maxbins[1] = 0;
		var ratio = 1;
		var maxwidth_set = 120;
		if(!outScopeflag){
			for(var i = 0 ; i < variants.length ; i++){
				if(variants[i].from > current_start){
					var natual_pos = ((variants[i].to+variants[i].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
					if(variants[i].genotypes[0] != undefined && variants[i].genotypes[0] != "0"){
						if(bins[0][Math.floor(natual_pos)] == undefined){
							bins[0][Math.floor(natual_pos)] = 1;
						}else{
							bins[0][Math.floor(natual_pos)] ++;
						}
						if(maxbins[0] < bins[0][Math.floor(natual_pos)]){
							maxbins[0] = bins[0][Math.floor(natual_pos)];
						}
					}
					if(variants[i].genotypes[1] != undefined && variants[i].genotypes[1] != "0"){
						if(bins[1][Math.floor(natual_pos)] == undefined){
							bins[1][Math.floor(natual_pos)] = 1;
						}else{
							bins[1][Math.floor(natual_pos)] ++;
						}
						if(maxbins[1] < bins[1][Math.floor(natual_pos)]){
							maxbins[1] = bins[1][Math.floor(natual_pos)];
						}
					}
				}
			}
		}else{
			var vlNodes = req3.responseXML.getElementsByTagName(xmlTagValues);
			for(var i=0; i<vlNodes.length; i++){
				if(vlNodes[i].getAttribute(xmlAttributeId).charAt(vlNodes[i].getAttribute(xmlAttributeId).length-1)=="L"){
					//alert("L");
					bins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
					maxbins[0] = Math.max.apply(null, bins[0]);
				}else{
					bins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
					maxbins[1] = Math.max.apply(null, bins[1]);
				}
			}
		}
		if(maxbins[0] > maxbins[1]){
			while((maxbins[0]*pixel_per_variant) > maxwidth_set){
				pixel_per_variant = pixel_per_variant/(pixel_per_variant + 1);
			}
		}else{
			while((maxbins[1]*pixel_per_variant) > maxwidth_set){
				pixel_per_variant = pixel_per_variant/(pixel_per_variant + 1);
			}
		}
		for(var i = 0 ; i < bins[0].length ; i++){
			if(bins[0][i] == undefined){
				bins[0][i] = 0;
			}
			R.rect(R_left+w/2-40-bins[0][i]*pixel_per_variant,i*bin_size+R_top,bins[0][i]*pixel_per_variant,bin_size,0).attr({fill:colO_hover2,stroke:colO_hover});
		}
		for(var i = 0 ; i < bins[1].length ; i++){
			if(bins[1][i] == undefined){
				bins[1][i] = 0;
			}
			R.rect(R_left+w/2+40,i*bin_size+R_top,bins[1][i]*pixel_per_variant,bin_size,0).attr({fill:colO_hover2,stroke:colO_hover});
		}
		var brwplotCanvas = $(document.getElementById("brwplot"));
		var markline = R.path("M"+R_left+","+R_top+" L"+R_left+","+(R_height-R_bottom)).attr({"stroke-dasharray": "-",stroke:"#000","stroke-width":1});
		var marklable = R.text(R_left+150, R_top-8, "").attr({font:"8px  \"Trebuchet MS\", Arial, sans-serif"});
		markline.hide();
		marklable.hide();
		var toucharea=[];
		toucharea[0] = R.rect(R_left+w/2-40-maxbins[0]*pixel_per_variant-6,R_top,maxbins[0]*pixel_per_variant+5,R_height-R_bottom-R_top,0).attr({stroke: "#000", fill: "#000", opacity: 0, cursor: "move"});
		toucharea[1] = R.rect(R_left+w/2+42,R_top,maxbins[1]*pixel_per_variant+6,R_height-R_bottom-R_top,0).attr({stroke: "#000", fill: "#000", opacity: 0, cursor: "move"});
		toucharea[0].mouseover(function(){
			var touch_x = restrictX - brwplotCanvas.position().left;
			var touch_y = restrictY - brwplotCanvas.position().top;
			var refSeqIndex = Math.round((R_left+w/2-40-touch_x)/pixel_per_variant);
			markline.stop().animate({transform:"t"+(touch_x-R_left)+" "+0},0);
			marklable.attr({text:refSeqIndex});
			markline.show();
			marklable.show();
		});
		toucharea[0].mousemove(function(){
			var touch_x = restrictX - brwplotCanvas.position().left;
			var touch_y = restrictY - brwplotCanvas.position().top;
			var refSeqIndex = Math.round((R_left+w/2-40-touch_x)/pixel_per_variant);
			markline.stop().animate({transform:"t"+(touch_x-R_left)+" "+0},0);
			marklable.attr({text:refSeqIndex});
			markline.show();
			marklable.show();
		});
		toucharea[0].mouseout(function(){
			markline.hide();
			marklable.hide();
		});
		toucharea[1].mouseover(function(){
			var touch_x = restrictX - brwplotCanvas.position().left;
			var touch_y = restrictY - brwplotCanvas.position().top;
			var refSeqIndex = Math.round((touch_x-R_left-w/2-40)/pixel_per_variant);
			markline.stop().animate({transform:"t"+(touch_x-R_left)+" "+0},0);
			marklable.attr({text:refSeqIndex});
			markline.show();
			marklable.show();
		});
		toucharea[1].mousemove(function(){
			var touch_x = restrictX - brwplotCanvas.position().left;
			var touch_y = restrictY - brwplotCanvas.position().top;
			var refSeqIndex = Math.round((touch_x-R_left-w/2-40)/pixel_per_variant);
			markline.stop().animate({transform:"t"+(touch_x-R_left)+" "+0},0);
			marklable.attr({text:refSeqIndex});
			markline.show();
			marklable.show();
		});
		toucharea[1].mouseout(function(){
			markline.hide();
			marklable.hide();
		});
	}else if(variants.length < l/font_height){
		lds_btnset.show();
		var current_pos = 0 - font_height;
		for(var i = 0 ; i < variants.length ; i++){
			if(variants[i].from > current_start){
				var natual_pos = ((variants[i].to+variants[i].from)/2 - current_start + 0.5)/(current_end - current_start + 1)*l;
				var name_pos = natual_pos;
				natual_pos += R_top;

				if(name_pos < current_pos + font_height){
					name_pos = current_pos + font_height;
				} else if(name_pos < i*font_height){
					name_pos = i*font_height;
				} else if(name_pos > l - (variants.length - i - 1)*font_height){
					name_pos = l - (variants.length - i - 1)*font_height;
				}	
				current_pos = name_pos;
				name_pos += R_top;
	
				variants[i].point = [];
				variants[i].name = [];
				variants[i].lines = [];

				var name_text = variants[i].id;

				if(variants[i].id.indexOf("Unnamed") == 0  && variants[i].dd && variants[i].dd.indexOf("rs") == 0){
					name_text = variants[i].dd + "*";
				}

				if(variants[i].genotypes[0] == undefined || variants[i].genotypes[0] == "0"){
					variants[i].point[0] = R.ellipse(R_left+w/2-30,natual_pos,2,2).attr({fill:"#FFF",stroke:colO});
				}else{
					variants[i].point[0] = R.ellipse(R_left+w/2-30,natual_pos,2,2).attr({fill:colO,stroke:colO});
					variants[i].lines[0] = R.path("M"+(R_left+text_length+2)+","+name_pos+" L"+(R_left+text_length+20)+","+name_pos+" L"+(R_left+w/2-52)+","+natual_pos+" L"+(R_left+w/2-33)+","+natual_pos).attr({stroke:colO,opacity:1});

					variants[i].name[0] = R.text(R_left+text_length,name_pos,name_text).attr({font:font_size2_text,opacity:1,"text-anchor":"end"}).attr({fill:colO,"font-weight":"normal"});

					(function(idx){
						variants[idx].name[0][0].style.cursor = "pointer";
						variants[idx].name[0][0].onmouseover = function(){
							variants[idx].name[0].animate({fill:colO_hover},100);
							var tran_x = variants[idx].point[0].attrs.cx - R_left;
							var tran_y = variants[idx].point[0].attrs.cy - R_top;
							floatbox[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
							floatbox[0].show();
							if(variants[idx].dd != undefined && variants[idx].dd != null){
								vdetail[0].attr({text:"dbSNPID: "+variants[idx].dd});
								vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
								vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
								if(variants[idx].type=="SNV"){
									vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
								}else{
									vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
								}
								vdetail[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
								floatboxset.stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
								vdetail[0].show();
								floatboxset.show();
							}else{
								vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
								vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
								if(variants[idx].type=="SNV"){
									vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
								}else{
									vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
								}
								floatboxset.stop().animate({transform:"t"+tran_x+" "+(tran_y-10)},0).toFront();
								floatboxset.show();
							}
							R.safari();
						};
						variants[idx].name[0][0].onmouseout = function(){
							variants[idx].name[0].animate({fill:colO},100);
							floatbox[0].hide();
							floatboxset.hide();
							vdetail[0].hide();
							R.safari();
						};
						variants[idx].name[0][0].onclick = function(){
							rightORleft = 0;
							select_a_variant(idx);
							R.safari();
						};
					})(i);
				}
				if(variants[i].genotypes[1] == undefined || variants[i].genotypes[1] == "0"){
					variants[i].point[1] = R.ellipse(R_left+w/2+30,natual_pos,2,2).attr({fill:"#FFF",stroke:colO});
				}else{
					variants[i].point[1] = R.ellipse(R_left+w/2+30,natual_pos,2,2).attr({fill:colO,stroke:colO});
					variants[i].lines[1] = R.path("M"+(R_left+w-text_length-2)+","+name_pos+" L"+(R_left+w-text_length-20)+","+name_pos+" L"+(R_left+w/2+52)+","+natual_pos+" L"+(R_left+w/2+33)+","+natual_pos).attr({stroke:colO,opacity:1});

					variants[i].name[1] = R.text(R_left+w-text_length,name_pos,name_text).attr({font:font_size2_text,opacity:1,"text-anchor":"start"}).attr({fill:colO});

					(function(idx){
						variants[idx].name[1][0].style.cursor = "pointer";
						variants[idx].name[1][0].onmouseover = function(){
							variants[idx].name[1].animate({fill:colO_hover},100);
							var tran_x = variants[idx].point[1].attrs.cx - R_left - 120 - 35;
							var tran_y = variants[idx].point[1].attrs.cy - R_top;
							floatbox[1].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
							floatbox[1].show();
							if(variants[idx].dd != undefined && variants[idx].dd != null){
								vdetail[0].attr({text:"dbSNPID: "+variants[idx].dd});
								vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
								vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
								if(variants[idx].type=="SNV"){
									vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
								}else{
									vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
								}
								vdetail[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
								floatboxset.stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
								vdetail[0].show();
								floatboxset.show();
							}else{
								vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
								vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
								if(variants[idx].type=="SNV"){
									vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
								}else{
									vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
								}
								floatboxset.stop().animate({transform:"t"+tran_x+" "+(tran_y-10)},0).toFront();
								floatboxset.show();
							}
							R.safari();
						};
						variants[idx].name[1][0].onmouseout = function(){
							variants[idx].name[1].animate({fill:colO},100);
							floatbox[1].hide();
							floatboxset.hide();
							vdetail[0].hide();
							R.safari();
						};
						variants[idx].name[1][0].onclick = function(){
							rightORleft = 1;
							select_a_variant(idx);
							R.safari();
						};
					})(i);
				}
			}
		}
		var fbox_w = 135;
		var fbox_h = 80;
		var fbox_arrowh = 10;
	
		floatbox[0] = R.path("M"+R_left+","+R_top+" L"+(R_left+fbox_arrowh)+","+(R_top-fbox_arrowh)+
					" L"+(R_left+fbox_arrowh)+","+(R_top-fbox_h/2)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_h/2)+
					" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_h/2)+" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_h/2)+
					" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_arrowh)+"z"
					).attr({stroke:colO,fill:"#ddd",stroke:colO_hover,"stroke-width":1,opacity:1});
		floatbox[1] = R.path("M"+(R_left+fbox_arrowh)+","+(R_top-fbox_h/2)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_h/2)+
					" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_arrowh)+" L"+(R_left+fbox_w+2*fbox_arrowh)+","+R_top+
				" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_arrowh)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_h/2)+
					" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_h/2)+"z"
					).attr({stroke:colO,fill:"#ddd",stroke:colO_hover,"stroke-width":1,opacity:1});
		floatbox[0].hide();
		floatbox[1].hide();
	
		var row_height = 15;
		var text_left = 15;
		var text_top = 15;
		var font_text_vdetail = "12px  \"Trebuchet MS\", Arial, sans-serif";
		vdetail[0] = R.text(R_left+text_left, R_top-fbox_h/2+text_top, "dbSNPID: ").attr({font:font_text_vdetail,"text-anchor":"start"});
		var fbox_w = 135;
		var fbox_h = 80;
		var fbox_arrowh = 10;
	
		floatbox[0] = R.path("M"+R_left+","+R_top+" L"+(R_left+fbox_arrowh)+","+(R_top-fbox_arrowh)+
					" L"+(R_left+fbox_arrowh)+","+(R_top-fbox_h/2)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_h/2)+
					" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_h/2)+" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_h/2)+
					" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_arrowh)+"z"
					).attr({stroke:colO,fill:"#ddd",stroke:colO_hover,"stroke-width":1,opacity:1});
		floatbox[1] = R.path("M"+(R_left+fbox_arrowh)+","+(R_top-fbox_h/2)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_h/2)+
					" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top-fbox_arrowh)+" L"+(R_left+fbox_w+2*fbox_arrowh)+","+R_top+
				" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_arrowh)+" L"+(R_left+fbox_arrowh+fbox_w)+","+(R_top+fbox_h/2)+
					" L"+(R_left+fbox_arrowh)+","+(R_top+fbox_h/2)+"z"
					).attr({stroke:colO,fill:"#ddd",stroke:colO_hover,"stroke-width":1,opacity:1});
		floatbox[0].hide();
		floatbox[1].hide();
		vdetail[0].hide();
		var row_height = 15;
		var text_left = 15;
		var text_top = 15;
		var font_text_vdetail = "12px  \"Trebuchet MS\", Arial, sans-serif";
		vdetail[0] = R.text(R_left+text_left, R_top-fbox_h/2+text_top, "dbSNPID: ").attr({font:font_text_vdetail,"text-anchor":"start"});
		vdetail[1] = R.text(R_left+text_left, R_top-fbox_h/2+text_top+row_height, "Pos: ").attr({font:font_text_vdetail,"text-anchor":"start"});
		vdetail[2] = R.text(R_left+text_left, R_top-fbox_h/2+text_top+2*row_height, "Genotype: ").attr({font:font_text_vdetail,"text-anchor":"start"});
		vdetail[3] = R.text(R_left+text_left, R_top-fbox_h/2+text_top+3*row_height, "Alt: ").attr({font:font_text_vdetail,"text-anchor":"start"});
	
		floatboxset = R.set();
		floatboxset.push(vdetail[1],vdetail[2],vdetail[3]);
		floatboxset.hide();
		vdetail[0].hide();
	} else if (variants.length < l/3){
		lds_btnset.hide();
		for(var i = 0 ; i < variants.length ; i++){
			if(variants[i].from > current_start){
				var natual_pos = ((variants[i].to+variants[i].from)/2 - current_start)/(current_end - current_start + 1)*l;
				natual_pos += R_top;

				variants[i].point = [];

				if(variants[i].genotypes[0] == undefined || variants[i].genotypes[0] == "0"){
					variants[i].point[0] = R.ellipse(R_left+w/2-30,natual_pos,2,2).attr({fill:"#FFF",stroke:colO});
				}else{
					variants[i].point[0] = R.ellipse(R_left+w/2-30,natual_pos,2,2).attr({fill:colO,stroke:colO});
				}
				if(variants[i].genotypes[1] == undefined || variants[i].genotypes[1] == "0"){
					variants[i].point[1] = R.ellipse(R_left+w/2+30,natual_pos,2,2).attr({fill:"#FFF",stroke:colO});
				}else{
					variants[i].point[1] = R.ellipse(R_left+w/2+30,natual_pos,2,2).attr({fill:colO,stroke:colO});
				}
			}
		}
	} 
}

function vcfTag2meaning(tag_value){
	var tag = "";
	var value = "";
	if(tag_value.indexOf("=")>0){
		tag = tag_value.split("=")[0];
		value = tag_value.split("=")[1];
	} else if(tag_value.indexOf(":")>0){
		tag = tag_value.split(":")[0];
		value = tag_value.split(":")[1];
	} else {
		tag = tag_value;
		value = null;
	}
	if(tag == "QUAL"){
		tag = "Quality";
	} else if(tag == "FILTER"){
		tag = "Filter";
	} else if(tag == "AA"){
		tag = "Ancestral Allele";
	} else if(tag == "AC"){
		tag = "Allele Count";
	} else if(tag == "AF"){
		tag = "Allele Frequency";
	} else if(tag == "AN"){
		tag = "Allele Number";
	} else if(tag == "BQ"){
		tag = "Base Quality";
	} else if(tag == "CIGAR"){
		tag = "CIGAR";
	} else if(tag == "H2"){
		tag = "HapMap2";
	} else if(tag == "H3"){
		tag = "HapMap3";
	} else if(tag == "MQ"){
		tag = "Mapping Quality";
	} else if(tag == "SB"){
		tag = "Strand Bias";
	} else if(tag == "SOMATIC"){
		tag = "Somatic Mutation";
	} else if(tag == "GMAF"){
		tag = "GMAF";
	} else if(tag == "GENEINFO"){
		tag = "Located Gene";
	} else if(tag == "dbSNPBuildID"){
		tag = "First dbSNP Build";
	} else if(tag == "NSF"){
		tag = "Non-synonymous Frameshift Allele";
	} else if(tag == "NSM"){
		tag = "Non-synonymous Missense Allele";
	} else if(tag == "NSN"){
		tag = "Non-synonymous Nonsense Allele";
	} else if(tag == "SYN"){
		tag = "Synonymous Allele";
	} else if(tag == "U3"){
		tag = "In 3'UTR Location";
	} else if(tag == "U5"){
		tag = "In 5'UTR Location";
	} else if(tag == "ASS"){
		tag = "In Acceptor Splice Site";
	} else if(tag == "DSS"){
		tag = "In Donor Splice Site";
	} else if(tag == "INT"){
		tag = "In Intron";
	} else if(tag == "R3"){
		tag = "In 3'Gene Region";
	} else if(tag == "R5"){
		tag = "In 5'Gene Region";
	} else if(tag == "MUT"){
		tag = "Mutation (Low Frequency)";
	} else if(tag == "G5A"){
		tag = "> 5% MAF In All Population";
	} else if(tag == "G5"){
		tag = "> 5% MAF In 1+ Population";
	} else if(tag == "GT"){
		tag = "Genotype";
	} else if(tag == "DP"){
		tag = "Read Depth";
	} else if(tag == "DS"){
		tag = "Genotype Dosage from MaCH/Thunder";
	} else if(tag == "FT"){
		tag = "Sample Genotype Filter";
	} else if(tag == "GL"){
		tag = "Genotype Likelihoods";
	} else if(tag == "GLE"){
		tag = "Genotype Likelihoods of Heterogenerous Ploidy";
	} else if(tag == "PL"){
		tag = "Phred-scaled Genotype Likelihoods";
	} else if(tag == "GP"){
		tag = "Phred-scaled Genotype Posterior Probabilities";
	} else if(tag == "GQ"){
		tag = "Conditional Genotype Quality";
	} else if(tag == "HQ"){
		tag = "Haplotype Qualities";
	} else if(tag == "PS"){
		tag = "Phase Set";
	} else if(tag == "PQ"){
		tag = "Phasing Quality";
	} else if(tag == "EC"){
		tag = "Expected Alternative Allele Counts";
	} else if(tag == "ID"){
		tag = "ID";
	} else if(tag == "Variant_seq"){
		tag = "Variant_seq";
	} else if(tag == "Genotype"){
		tag = "Genotype";
	} else if(tag == "Zygosity"){
		tag = "Zygosity";
	} else if(tag == "Dbxref"){
		tag = "Dbxref";
	} else if(tag == "Reference_seq"){
		tag = "Reference_seq";
	} else if(tag.indexOf("_AF")>0){
		tag = tag;
	} else {
		tag = null;
	}
	var infohtml = "";
	if(tag != null && value !=null){
		infohtml = "<tr><td>"+tag+":</td><td>"+value+"</td></tr>";
	} else if (tag != null){
		infohtml = "<tr><td colspan=\"2\">"+tag+"</td></tr>";
	}
	return infohtml;
}

function show_vardetail(){
	if(reqVD.readyState == 4) {
		if(reqVD.status == 200){
			var XMLDoc = reqVD.responseXML;
			var variantNode =  XMLDoc.getElementsByTagName(xmlTagVariant)[0];
			var variantId = variantNode.getAttribute(xmlAttributeId);
			var variantType = variantNode.getAttribute(xmlAttributeType);
			var variantFrom = variantNode.getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
			var variantTo = variantNode.getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
			var variantLetter = null;
			var variant_dbSNPID = null;
			var variant_BLS_mate = null;
			if(variantNode.getElementsByTagName(xmlTagLetter).length > 0){
				variantLetter = variantNode.getElementsByTagName(xmlTagLetter)[0].childNodes[0].nodeValue;
			}
			if(variantNode.getAttribute(xmlAttribute_dbSNPID)){
				variant_dbSNPID = variantNode.getAttribute(xmlAttribute_dbSNPID);
			}
			if(variantNode.getElementsByTagName(xmlTagMate).length > 0){
				variant_BLS_mate = variantNode.getElementsByTagName(xmlTagMate)[0].childNodes[0].nodeValue;
			}
			var variantDescription = variantNode.getElementsByTagName(xmlTagDescription)[0].childNodes[0].nodeValue + "";

			var vdtable=document.getElementById("vdtablebody");
			var varname = variants[csv].id;
			if(variants[csv].id.indexOf("Unnamed") == 0  && variants[csv].dd && variants[csv].dd.indexOf("rs") == 0){
				varname = variants[csv].dd + "*";
			}
			document.getElementById("varname").innerHTML = varname;
			if(variantId!="."&&(/^rs/).test(variantId)){
				$("#vardetail_link").css("display","block");
				document.getElementById("vardetail_link").onclick = function(){
					window.open(refSNPurl+variantId);
				};
			}else if(variant_dbSNPID&&(/^rs/).test(variant_dbSNPID)){
				$("#vardetail_link").css("display","block");
				document.getElementById("vardetail_link").onclick = function(){
					window.open(refSNPurl+variant_dbSNPID);
				};
			}else{
				$("#vardetail_link").css("display","none");
				//document.getElementById("vardetail_link").innerHTML = "";
				document.getElementById("vardetail_link").onclick = function(){};
			}
			document.getElementById("vardetail_link").onmouseover = function(){
				document.getElementById("vardetail_link").style.color="#FF0";
			};
			document.getElementById("vardetail_link").onmouseout = function(){
				document.getElementById("vardetail_link").style.color="blue";
			};
			document.getElementById("vardetail_scale").innerHTML = variants[csv].chr+":" + variantFrom + "-" + variantTo;
			document.getElementById("vardetail_type").innerHTML = variantType;

			if(variantLetter){
				$("#vardetail_letter_trNode").css("display","table-row");
				$("#vardetail_letter").html(variantLetter);
			}else{
				$("#vardetail_letter_trNode").css("display","none");
			}
			if(variant_dbSNPID){
				$("#vardetail_dbSNPID_trNode").css("display","table-row");
				$("#vardetail_dbSNPID").html(variant_dbSNPID);
			}else{
				$("#vardetail_dbSNPID_trNode").css("display","none");
			}
			if(variant_BLS_mate){
				$("#vardetail_Mate_trNode").css("display","table-row");
				$("#vardetail_Mate").html(variant_BLS_mate);
			}else{
				$("#vardetail_Mate_trNode").css("display","none");
			}

			var infohtml = "";
			var variantInfos = variantDescription.split(";");
			infohtml+="<tr bgcolor=\"#ccc\"><td colspan=\"2\" width=\"193\"><b>Detail in VCF fields:</b></td></tr>";
			for(var infoidx = 0;infoidx < variantInfos.length;infoidx++){
				infohtml+=vcfTag2meaning(variantInfos[infoidx]);
			}
			if(variantNode.getElementsByTagName("SampleInfo").length > 0){
				infohtml+="<tr bgcolor=\"#ccc\"><td colspan=\"2\"><b>Detail in Sample field:</b></td></tr>";
				var variantFormat = variantNode.getElementsByTagName("SampleInfo")[0].childNodes[0].nodeValue + "";
				var variantSinfos= variantFormat.split(";");
				for(var infoidx = 0;infoidx < variantSinfos.length;infoidx++){
						infohtml+=vcfTag2meaning(variantSinfos[infoidx]);
				}
			}
			if(variantNode.getElementsByTagName("Dbsnp").length > 0){
				infohtml+="<tr bgcolor=\"#ccc\"><td colspan=\"2\"><b>Detail in dbSNP:</b></td></tr>";
				var variantDbsnp = variantNode.getElementsByTagName("Dbsnp")[0].childNodes[0].nodeValue + "";
				var variantDinfos= variantDbsnp.split(";");
				for(var infoidx = 0;infoidx < variantDinfos.length;infoidx++){
						infohtml+=vcfTag2meaning(variantDinfos[infoidx]);
				}
			}

			$("#vardetail_info_table").html(infohtml);

			var brwplotCanvas = $(document.getElementById("brwplot"));
			var tablewidth = 220;
			var tableheight = 330;
			$(document.getElementById("vardetail")).css("position", "absolute");
			
			if(variants[csv].name[rightORleft].attrs.x < R_width - tablewidth){
				$(document.getElementById("vardetail")).css("left", brwplotCanvas.position().left + variants[csv].name[rightORleft].attrs.x + rightORleft*50 - 64);
			}else{
				$(document.getElementById("vardetail")).css("left", brwplotCanvas.position().left + variants[csv].name[rightORleft].attrs.x - tablewidth + rightORleft*50);
			}
			if(variants[csv].name[rightORleft].attrs.y + 5 < R_height - tableheight && brwplotCanvas.position().top + variants[csv].name[rightORleft].attrs.y + tableheight < document.body.clientHeight+document.body.scrollTop){
				$(document.getElementById("vardetail")).css("top", brwplotCanvas.position().top + variants[csv].name[rightORleft].attrs.y + 5);
			}else{
				$(document.getElementById("vardetail")).css("top", brwplotCanvas.position().top + variants[csv].name[rightORleft].attrs.y - 16 - tableheight);
			}
			
			/*if(variants[csv].name[rightORleft][0].offsetLeft < R_width - tablewidth){
					$(document.getElementById("vardetail")).css("left", brwplotCanvas.position().left + variants[csv].name[rightORleft][0].offsetLeft + rightORleft*50);
			}else{
					$(document.getElementById("vardetail")).css("left", brwplotCanvas.position().left + variants[csv].name[rightORleft][0].offsetLeft - tablewidth + rightORleft*50);
			}
			if(variants[csv].name[rightORleft][0].offsetTop + 15 < R_height - tableheight && brwplotCanvas.position().top + variants[csv].name[rightORleft][0].offsetTop + tableheight < document.body.clientHeight+document.body.scrollTop){
					$(document.getElementById("vardetail")).css("top", brwplotCanvas.position().top + variants[csv].name[rightORleft][0].offsetTop + 15);
			}else{
					$(document.getElementById("vardetail")).css("top", brwplotCanvas.position().top + variants[csv].name[rightORleft][0].offsetTop - 10 - tableheight);
			}*/
			
			if(document.getElementById("vardetail").style.display == "none"){
				$(document.getElementById("vardetail")).css("display", "block");
			}else{
				close_vardetail();
			}
		}
	}
}

function show_same_collectionvar(){
	if(reqV.readyState == 4) {
		if(reqV.status == 200){
			var xmlDoc=null;
			var xmlString = reqV.responseText;
			if(!window.DOMParser && window.ActiveXObject){
				var xmlDomVersions = ['MSXML.2.DOMDocument.6.0','MSXML.2.DOMDocument.3.0','Microsoft.XMLDOM'];
				for(var i=0;i<xmlDomVersions.length;i++){
					try{
						xmlDoc = new ActiveXObject(xmlDomVersions[i]);
						xmlDoc.async = false;
						xmlDoc.loadXML(xmlString); 
						break;
					}catch(e){
					}
				}
			}
			else if(window.DOMParser && document.implementation && document.implementation.createDocument){
				try{
					domParser = new  DOMParser();
					xmlDoc = domParser.parseFromString(xmlString, 'text/xml');
				}catch(e){
				}
			}
			for(var i in individuals){
				individuals[i].gtmark = "0|0";
			}
			if(csi!=undefined){	
				for(var family_member in families){
					if(family_member != "roots" && individuals[family_member].ifs == "true" && families[family_member].markobj != undefined){
						for(var temp_root in families[family_member].markobj){
							//families[family_member].markobj[temp_root].hide();
							families[family_member].markobj[temp_root].attr({text:"0|0",fill:"#000"});
						}
					}
				}
				if(xmlDoc.getElementsByTagName("IndividualOrder")[0]!=null && xmlDoc.getElementsByTagName("IndividualOrder")[0].firstChild!=null){			
				//	alert("OK, the value is " + xmlDoc.getElementsByTagName("IndividualOrder")[0].firstChild.nodeValue);
					var collectionlist = xmlDoc.getElementsByTagName("IndividualOrder")[0].firstChild.nodeValue.split(",");
					if(collectionlist!=null){
						for(var j in collectionlist){
							var id_info = collectionlist[j].split(":");
							if(individuals[id_info[0]] != undefined){
								individuals[id_info[0]].gtmark = id_info[2];
							}
							if(individuals[id_info[0]].ifs == "true"){
								for(var temp_root in families[id_info[0]].markobj){
									if(families[id_info[0]].markobj[temp_root] != undefined){
										families[id_info[0]].markobj[temp_root].attr({text:id_info[2],fill:"#F00"});
										families[id_info[0]].markobj[temp_root].show();
									}
								}														
							}
						}
					}
				}
				else{
//					alert("Nothing!");
				}
			}
			list_individuals();
		}
	}
}

function select_a_variant(idx){
	if(variants[idx] != undefined){
		////////////////////insert request of searchIndividual
		/*alert(variants[idx].chr + ":"  + variants[idx].from + ":" + variants[idx].to + ":" + variants[idx].letter);
		var varString;
		if(variants[idx].letter=="--"){
			varString = variants[idx].chr + ":"  + variants[idx].from + ":" + variants[idx].to + ":-";
		}else{
			varString = variants[idx].chr + ":"  + variants[idx].from + ":" + variants[idx].to + ":" + variants[idx].letter;
		}
		var formV = new FormData();
		formV.append("variants",varString);
		formV.append("enctype","multipart/form-data");
		reqV.onreadystatechange = show_same_collectionvar;
		querry = "action=searchIndividual&tracks="+trackname;
		reqV.open("POST","servlet/test.do?"+querry,true);
		reqV.send(formV);
		
		*//////////////////////////////

		var radioObj = document.getElementById(idx+"__varlist_select")
		if(variants[idx].selected){
			for(var family_member in families){
				if(family_member != "roots" && individuals[family_member].ifs == "true" && families[family_member].markobj != undefined){
					for(var temp_root in families[family_member].markobj){
						families[family_member].markobj[temp_root].attr({text:""});
					}
				}
			}

			for(var i in individuals){
				individuals[i].gtmark = "--";
			}
			variants[idx].selected = false;
			radioObj.checked = false;
			if(variants[idx].name != undefined){
				if(variants[idx].genotypes[0] != undefined && variants[idx].genotypes[0] != "0"){
					variants[idx].name[0].attr({"font-weight":"normal"});
					variants[idx].point[0].attr({"stroke-width":1});
					variants[idx].lines[0].attr({"stroke-width":1});
				}
				if(variants[idx].genotypes[1] != undefined && variants[idx].genotypes[1] != "0"){
					variants[idx].name[1].attr({"font-weight":"normal"});
					variants[idx].point[1].attr({"stroke-width":1});
					variants[idx].lines[1].attr({"stroke-width":1});
				}
				if(variants[idx].functional_v != undefined){
					for(var fvidx = 0 ; fvidx < variants[idx].functional_v.length ; fvidx++){
						variants[idx].functional_v[fvidx].obj.attr({"font-weight":"normal","stroke-width":1});
					}
				}
				if(variants[idx].lds != undefined){
					variants[idx].lds.hide();
				}
				if(lds_flag == 1){
					all_lds.show();
				}else{
					all_lds.hide();
				}
			}
			csv = -1;
		} else {
			var varString;
			if(variants[idx].letter=="--"){
				varString = variants[idx].chr + ":"  + variants[idx].from + ":" + variants[idx].to + ":-";
			}else{
				varString = variants[idx].chr + ":"  + variants[idx].from + ":" + variants[idx].to + ":" + variants[idx].letter;
			}
			var formV = new FormData();
			formV.append("variants",varString);
			formV.append("enctype","multipart/form-data");
			reqV.onreadystatechange = show_same_collectionvar;
			querry = "action=searchIndividual&tracks="+trackname;
			reqV.open("POST","servlet/test.do?"+querry,true);
			reqV.send(formV);

			if(csv>=0 && variants[csv] != undefined){
				variants[csv].selected = false;
			
				if(variants[csv].name != undefined){
					if(variants[csv].genotypes[0] != undefined && variants[csv].genotypes[0] != "0"){
						variants[csv].name[0].attr({"font-weight":"normal"});
						variants[csv].point[0].attr({"stroke-width":1});
						variants[csv].lines[0].attr({"stroke-width":1});
					}
					if(variants[csv].genotypes[1] != undefined && variants[csv].genotypes[1] != "0"){
						variants[csv].name[1].attr({"font-weight":"normal"});
						variants[csv].point[1].attr({"stroke-width":1});
						variants[csv].lines[1].attr({"stroke-width":1});
					}
					if(variants[csv].functional_v != undefined){
						for(var fvidx = 0 ; fvidx < variants[csv].functional_v.length ; fvidx++){
							variants[csv].functional_v[fvidx].obj.attr({"font-weight":"normal","stroke-width":1});
						}
					}
					if(variants[csv].lds != undefined){
						variants[csv].lds.hide();
					}
				}
			}
			variants[idx].selected = true;
			radioObj.checked = true;
			if(variants[idx].name != undefined){
				if(variants[idx].genotypes[0] != undefined && variants[idx].genotypes[0] != "0"){
					variants[idx].name[0].attr({"font-weight":"bolder"});
					variants[idx].point[0].attr({"stroke-width":2});
					variants[idx].lines[0].attr({"stroke-width":2});
				}
				if(variants[idx].genotypes[1] != undefined && variants[idx].genotypes[1] != "0"){
					variants[idx].name[1].attr({"font-weight":"bolder"});
					variants[idx].point[1].attr({"stroke-width":2});
					variants[idx].lines[1].attr({"stroke-width":2});
				}
				if(variants[idx].functional_v != undefined){
					for(var fvidx = 0 ; fvidx < variants[idx].functional_v.length ; fvidx++){
						variants[idx].functional_v[fvidx].obj.attr({"font-weight":"bolder","stroke-width":2});
					}
				}
				all_lds.hide();
				if(variants[idx].lds != undefined && lds_flag == 1){
					variants[idx].lds.show();
				}
			}
			csv = idx;
			if(document.getElementById("varlist").style.display == "none"){
				reqVD.onreadystatechange = show_vardetail;
				var querryVD;
				if(variants[idx].id.indexOf("Unnamed") == 0){
					querryVD = "action=getDetail&tracks=_"+csi+"@"+trackname+"&id=.&start="+variants[idx].from+"&end="+variants[idx].to;
				}else{
					querryVD = "action=getDetail&tracks=_"+csi+"@"+trackname+"&id="+variants[idx].id+"&start="+variants[idx].from+"&end="+variants[idx].to;
				}
				reqVD.open("GET","servlet/test.do?"+querryVD,false);
				reqVD.send();
				floatbox[0].hide();
				floatbox[1].hide();
				floatboxset.hide();
				vdetail[0].hide();
			}
		}
	}
	list_individuals();
}

function change_variant_color(vPointer,color){
	var l = R_height - R_top - R_bottom;
	var w = R_width - R_left - R_right;
	if(variants[vPointer].name != undefined){

		var color2 = colO_hover;
		if(color == colP){
			color2 = colP_hover;
		}else if(color == colD){
			color2 = colD_hover;
		}else if(color == colM){
			color2 = colM_hover;
		}else if(color == colInter){
			color2 = colInter_hover;
		}else if(color == colDiff){
			color2 = colDiff_hover;
		}

		if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
			variants[vPointer].name[0].attr({fill:color}).toFront();

			(function(idx,color){
				variants[idx].name[0][0].style.cursor = "pointer";
				variants[idx].name[0][0].onmouseover = function(){
					variants[idx].name[0].animate({fill:color2},100);
					
					var tran_x = variants[idx].point[0].attrs.cx - R_left;
					var tran_y = variants[idx].point[0].attrs.cy - R_top;
					floatbox[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
					floatbox[0].show();
					if(variants[idx].dd != undefined && variants[idx].dd != null){
						vdetail[0].attr({text:"dbSNPID: "+variants[idx].dd});
						vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
						vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
						if(variants[idx].type=="SNV"){
							vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
						}else{
							vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
						}
						vdetail[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
						floatboxset.stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
						vdetail[0].show();
						floatboxset.show();
					}else{
						vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
						vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
						if(variants[idx].type=="SNV"){
							vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
						}else{
							vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
						}
						floatboxset.stop().animate({transform:"t"+tran_x+" "+(tran_y-10)},0).toFront();
						floatboxset.show();
					}
					R.safari();
				};
				variants[idx].name[0][0].onmouseout = function(){
					variants[idx].name[0].animate({fill:color},100);
					floatbox[0].hide();
					floatboxset.hide();
					vdetail[0].hide();
					R.safari();
				};
				variants[idx].name[0][0].onclick = function(){
					rightORleft = 0;
					select_a_variant(idx);
					R.safari();
				};
			})(vPointer,color);
		}
		if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
			variants[vPointer].name[1].attr({fill:color}).toFront();

			(function(idx,color){
				variants[idx].name[1][0].style.cursor = "pointer";
				variants[idx].name[1][0].onmouseover = function(){
					variants[idx].name[1].animate({fill:color2},100);
					var tran_x = variants[idx].point[1].attrs.cx - R_left - 120 - 35;
					var tran_y = variants[idx].point[1].attrs.cy - R_top;
					floatbox[1].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
					floatbox[1].show();
					if(variants[idx].dd != undefined && variants[idx].dd != null){
						vdetail[0].attr({text:"dbSNPID: "+variants[idx].dd});
						vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
						vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
						if(variants[idx].type=="SNV"){
							vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
						}else{
							vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
						}
						vdetail[0].stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
						floatboxset.stop().animate({transform:"t"+tran_x+" "+tran_y},0).toFront();
						vdetail[0].show();
						floatboxset.show();
					}else{
						vdetail[1].attr({text:"Position: "+current_chr+":"+variants[idx].from});
						vdetail[2].attr({text:"Genotype: "+variants[idx].genotype});
						if(variants[idx].type=="SNV"){
							vdetail[3].attr({text:"Alt allele: "+variants[idx].letter});
						}else{
							vdetail[3].attr({text:"Alt allele: <"+variants[idx].type+">"});
						}
						floatboxset.stop().animate({transform:"t"+tran_x+" "+(tran_y-10)},0).toFront();
						floatboxset.show();
					}
					R.safari();
				};
				variants[idx].name[1][0].onmouseout = function(){
					variants[idx].name[1].animate({fill:color},100);
					floatbox[1].hide();
					floatboxset.hide();
					vdetail[0].hide();
					R.safari();
				};
				variants[idx].name[1][0].onclick = function(){
					rightORleft = 1;
					select_a_variant(idx);
					R.safari();
				};
			})(vPointer,color);
		}
	}
	if(variants[vPointer].point != undefined && variants.length < l/3){
		if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
			variants[vPointer].point[0].attr({fill:color,stroke:color}).toFront();
		}
		if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
			variants[vPointer].point[1].attr({fill:color,stroke:color}).toFront();
		}
	} else if(color != colO && variants.length < l/3){
		var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*l + R_top;
		variants[vPointer].point = [];
		if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
			variants[vPointer].point[0] = R.ellipse(R_left+w/2-30,natual_pos,2,2).attr({fill:color,stroke:color});
		}
		if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
			variants[vPointer].point[1] = R.ellipse(R_left+w/2+30,natual_pos,2,2).attr({fill:color,stroke:color});
		}
	} else if(variants.length >= l/3){
		if(variants[vPointer].point != undefined){
			if(variants[vPointer].point[0] != undefined){
				variants[vPointer].point[0].remove();
			}
			if(variants[vPointer].point[1] != undefined){
				variants[vPointer].point[1].remove();
			}
			delete variants[vPointer].point;
		}
	}
	if(variants[vPointer].lines != undefined){
		if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
			variants[vPointer].lines[0].attr({stroke:color}).toFront();
		}
		if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
			variants[vPointer].lines[1].attr({stroke:color}).toFront();
		}
	}
}

function close_gdflist(){
	$(document.getElementById("gdflist")).css("display", "none");
}

function close_vardetail(){
	$(document.getElementById("vardetail")).css("display", "none");
}

function plot_genes(){
	var l = R_height - R_top - R_bottom;
	var w = R_width - R_left - R_right;
	var font_height = 15;
	var box_width = 10;
	var band_width = 6;
	var deslinelength = 80;
	var brwplotCanvas = $(document.getElementById("brwplot"));	

	var vari = 0;
	var plotted_amino_var = {};
	var current_pos = 0 - font_height;
	for(var i = 0 ; i < functional_vPointer.length ; i++){
		for(var fvidx = 0 ; fvidx < variants[functional_vPointer[i]].functional_v.length ; fvidx++){
			var var_top = map_coord(l,variants[functional_vPointer[i]].functional_v[fvidx].from);
			var var_bot = map_coord(l,variants[functional_vPointer[i]].functional_v[fvidx].to);
			if(var_top == l || var_bot == 0){
				continue;
			}
			var natual_pos = (var_top+var_bot)/2;
			var name_pos = natual_pos;
			natual_pos += R_top;
			if(name_pos < current_pos + font_height){
				name_pos = current_pos + font_height;
			} else if(name_pos < vari*font_height){
				name_pos = vari*font_height;
			} else if(name_pos > l - (functional_vnum - vari - 1)*font_height){
				name_pos = l - (functional_vnum - vari - 1)*font_height;
			}
			vari++;

			current_pos = name_pos;
			name_pos += R_top;
			variants[functional_vPointer[i]].functional_v[fvidx].obj = draw_var(natual_pos,name_pos,l,box_width,variants[functional_vPointer[i]].functional_v[fvidx].letter);
		}
	}
	for(var i = 0 ; i < genes.length ; i++){
		genes[i].obj = R.set();
		
		if(genes[i].subs == undefined){
			var main = draw_box(genes[i].from,genes[i].to,l,box_width);
			genes[i].obj.push(main);
//			var direction = draw_triangle(genes[i].from,genes[i].to,l,genes[i].strand,box_width);
//			if(direction != null){
//				genes[i].obj.push(direction);
//			}
		} else {
			for(var j = 0 ; j < genes[i].subs.length ; j++){
				var temp = null;
				if(genes[i].subs[j].type == "X"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "shBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "eBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "skBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "lBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "psBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "seBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "pseBox"){
					temp = draw_box(genes[i].subs[j].from,genes[i].subs[j].to,l,box_width);
				} else if(genes[i].subs[j].type == "D"){
					temp = draw_band(genes[i].subs[j].from,genes[i].subs[j].to,l,band_width,box_width);
				} else if(genes[i].subs[j].type == "eBand"){
					temp = draw_band(genes[i].subs[j].from,genes[i].subs[j].to,l,band_width,box_width);
				} else if(genes[i].subs[j].type == "sBand"){
					temp = draw_band(genes[i].subs[j].from,genes[i].subs[j].to,l,band_width,box_width);
				} else if(genes[i].subs[j].type == "L"){
					temp = draw_line(genes[i].subs[j].from,genes[i].subs[j].to,l,genes[i].strand,box_width);
				}
				if(temp != null){
					genes[i].obj.push(temp);
				}
			}
		}
		if(genes[i].vars != undefined){
			for(var fv in genes[i].vars){
				genes[i].obj.push(variants[functional_v[fv]].functional_v[genes[i].vars[fv]].obj);
			}
		}
	}
	if(genes.length > 0){
		var font_size2_text = "14px \"Trebuchet MS\", Arial, sans-serif";
		draw_arrow(R_width-R_right+10+box_width/2,0,-1);
		draw_arrow(R_width-R_right+10+box_width/2,0,1);
		draw_arrow(R_width-R_right+10+box_width/2,1,-1);
		draw_arrow(R_width-R_right+10+box_width/2,1,1);
		cssObj = R.text(R_width-R_right+10+box_width/2,R_top-8-20,"ALL").attr({font:font_size2_text});
		cstObj = R.text(R_width-R_right+10+box_width/2,R_top-8,"ALL").attr({font:font_size2_text});
	}
	
	var w = 55;
	var h = 25;
	var xm = R_width + 70;
	var y = (R_height-R_bottom-R_top)/2;
	gdf_icon = R.set();
	var font_size2_text = "15px \"Trebuchet MS\", Arial, sans-serif";
	//gdf_icon = R.path("M26.711,14.086L16.914,4.29c-0.778-0.778-2.051-0.778-2.829,0L4.29,14.086c-0.778,0.778-0.778,2.05,0,2.829l9.796,9.796c0.778,0.777,2.051,0.777,2.829,0l9.797-9.797C27.488,16.136,27.488,14.864,26.711,14.086zM14.702,8.981c0.22-0.238,0.501-0.357,0.844-0.357s0.624,0.118,0.844,0.353c0.221,0.235,0.33,0.531,0.33,0.885c0,0.306-0.101,1.333-0.303,3.082c-0.201,1.749-0.379,3.439-0.531,5.072H15.17c-0.135-1.633-0.301-3.323-0.5-5.072c-0.198-1.749-0.298-2.776-0.298-3.082C14.372,9.513,14.482,9.22,14.702,8.981zM16.431,21.799c-0.247,0.241-0.542,0.362-0.885,0.362s-0.638-0.121-0.885-0.362c-0.248-0.241-0.372-0.533-0.372-0.876s0.124-0.638,0.372-0.885c0.247-0.248,0.542-0.372,0.885-0.372s0.638,0.124,0.885,0.372c0.248,0.247,0.372,0.542,0.372,0.885S16.679,21.558,16.431,21.799z").attr({fill: colO_hover, stroke: "none",cursor:"pointer"});	
	//gdf_icon.stop().animate({transform: "t"+(xm-65)+" "+y+" s2.0"});	
	//gdf_icon.push(R.rect(xm-27,y-12,w,h,2).attr({fill:colO_hover,stroke:colO_hover}));
	//var temp_gdficon = R.rect(xm-22,y-7,w-10,h-10,2).attr({fill:"#FFF",stroke:colO_hover});	
	//gdf_icon.push(temp_gdficon);
	var temp_gdficon = R.rect(xm-30,y-10,w+4,h+10,3).attr({fill:colO_hover,cursor:"pointer",stroke:colO_hover,opacity:0.5});
	gdf_icon.push(R.text(xm,y,"OMIM").attr({font:font_size2_text,cursor:"pointer"}));
	gdf_icon.push(R.text(xm,y+15,"Records").attr({font:font_size2_text,cursor:"pointer"}));
	gdf_icon.push(temp_gdficon);
	
	gdf_icon.hide();
	if(gdf_list.length > 0){
		gdf_icon.show();
		gdf_icon.mouseover(function(){
			temp_gdficon.animate({fill:colO_hover2},200);
		});
		gdf_icon.mouseout(function(){
			temp_gdficon.animate({fill:colO_hover},200);
		});
		gdf_icon.mousedown(function(){
			var t=document.getElementById("gdftablebody");
			t.innerHTML = "";
			for(var i=0; i<(gdf_list.length+1); i++){
				var row=document.createElement("tr");
				row.id = "row"+i;
				var cell3=document.createElement("td");
				var cell2=document.createElement("td");
				var cell1=document.createElement("td");
				if(i==0){
					cell3.appendChild(document.createTextNode("Source"));
					cell2.appendChild(document.createTextNode("Symbol"));
					cell1.appendChild(document.createTextNode("Id"));
				}else{
					cell3.id="gdfsource"+i;
					cell3.appendChild(document.createTextNode(gdf_list[i-1].source));
					cell2.appendChild(document.createTextNode(gdf_list[i-1].symbol));
					cell1.appendChild(document.createTextNode(gdf_list[i-1].id));
				}
				row.appendChild(cell1);
				row.appendChild(cell2);
				row.appendChild(cell3);
				document.getElementById("gdftablebody").appendChild(row);
			}
			var rows=document.getElementById("gdf_list").rows;  
			document.getElementById("gdf_list").style.fontSize="13px";  
			if(rows.length>1){  
				for(var i=1;i<rows.length;i++){  
					(function(i){  
						var obj=rows[i];
						var patternOMIM = /^.*,.*(\d{6}).*$/;
						if(patternOMIM.test(gdf_list[i-1].id)){
							$(document.getElementById("gdfsource"+i)).css("cursor","pointer"); 
							document.getElementById("gdfsource"+i).style.textDecoration="underline";
							var OMIMentry = RegExp.$1;
							obj.cells[2].onclick = function(){
								window.open(OMIMurl+OMIMentry);
							};
							obj.cells[2].onmouseover = function(){
								obj.cells[2].style.color="#FF0";
							};
							obj.cells[2].onmouseout = function(){
								obj.cells[2].style.color="#FFF";
							};
						}
					})(i);  
				}  
			} 
			$(document.getElementById("gdflist")).css("position", "absolute");
			$(document.getElementById("gdflist")).css("top", brwplotCanvas.position().top+y+27);
			$(document.getElementById("gdflist")).css("left", brwplotCanvas.position().left+xm-460);
			if(document.getElementById("gdflist").style.display == "none"){
				$(document.getElementById("gdflist")).css("display", "block");
			}else{
				close_gdflist();
			}
		});
	}
	document.body.addEventListener("mousedown", mousedownOutsideRepeatTooltip, false);
	function mousedownOutsideRepeatTooltip(evt){
		evt = evt || window.event;
		var eventTarget = evt.target || evt.srcElement;
		var flagomim=0;
		var flagvar=0;
		while(eventTarget){
			if(eventTarget==document.getElementById("gdflist")|| eventTarget==gdf_icon[0][0] || eventTarget==gdf_icon[1][0] || eventTarget==gdf_icon[2][0]){
				flagomim=1;
			}
			if(eventTarget==document.getElementById("vardetail")){
				flagvar=1;
			}
			eventTarget = eventTarget.parentNode;
		}
		if(flagomim==0) {
			document.getElementById("gdflist").style.display = "none";
		}
		if(flagvar==0) {
			document.getElementById("vardetail").style.display = "none";
		}
	}
	function change_symbol(x,sh,direction){
		var font_size2_text = "14px \"Trebuchet MS\", Arial, sans-serif";
		if(cssObj == null){
			return;
		}
		cssObj.remove();
		cstObj.remove();
		if(css == 0 && direction == -1){
			css = symbols.length - 1;
		} else if(css == symbols.length -1 && direction == 1){
			css = 0;
		} else {
			css += direction;
		}
		if(css == 0){
			for(var i = 1 ; i < symbols.length ; i++){
				for(var j = 1 ; j < symbols[i].length ; j++){
					genes[symbols[i][j]].obj.show();
				}
			}
			cssObj = R.text(x,R_top-sh-20,"ALL").attr({font:font_size2_text});
		} else {
			for(var i = 1 ; i < symbols.length ; i++){
				for(var j = 1 ; j < symbols[i].length ; j++){
					genes[symbols[i][j]].obj.hide();
				}
			}
			if(css > 0 && css < symbols.length){
				for(var j = 1 ; j < symbols[css].length ; j++){
					genes[symbols[css][j]].obj.show();
				}
			}
			cssObj = R.text(x,R_top-sh-20,symbols[css][0]).attr({font:font_size2_text});
			cssObj[0].style.cursor = "pointer";
			cssObj[0].onclick = function (){
				var reqHGNC = createXMLHttpRequest();
				reqHGNC.open("GET","servlet/test.do?action=getGene&gene="+symbols[css][0],false);
				reqHGNC.send(null);
				var hgnc_id = reqHGNC.responseXML.getElementsByTagName("HGNC")[0].childNodes[0].nodeValue;
				window.open(hgnc_url+hgnc_id);
			};
			cssObj[0].onmouseover = function (){
				cssObj[0].style.textDecoration = "underline";
				cssObj.attr({fill:colO_hover2});
			};
			cssObj[0].onmouseout = function (){
			//	cssObj.attr({"font-weight":"normal"});
				cssObj[0].style.textDecoration = "";
				cssObj.attr({fill:"#000"});
			};
		}

		cst = 0;
		cstObj = R.text(x,R_top-sh,"ALL").attr({font:font_size2_text});
	}
	function change_transcript(x,sh,direction){
		var font_size2_text = "14px \"Trebuchet MS\", Arial, sans-serif";
		var font_size1_text = "13px \"Trebuchet MS\", Arial, sans-serif";
		if(cssObj == null || cstObj == null || css == 0){
			return;
		}
		cstObj.remove();
		if(cst == 0 && direction == -1){
			cst = symbols[css].length - 1;
		} else if(cst == symbols[css].length -1 && direction == 1){
			cst = 0;
		} else {
			cst += direction;
		}
		if(cst == 0){
			for(var i = 1 ; i < symbols.length ; i++){
				for(var j = 1 ; j < symbols[i].length ; j++){
					genes[symbols[i][j]].obj.hide();
				}
			}
			if(css > 0 && css < symbols.length){
				for(var j = 1 ; j < symbols[css].length ; j++){
					genes[symbols[css][j]].obj.show();
				}
			}
			cstObj = R.text(x,R_top-sh,"ALL").attr({font:font_size2_text});
		} else {
			for(var i = 1 ; i < symbols.length ; i++){
				for(var j = 1 ; j < symbols[i].length ; j++){
					genes[symbols[i][j]].obj.hide();
				}
			}
			if(css > 0 && css < symbols.length 
			&& cst > 0 && cst < symbols[css].length){
				genes[symbols[css][cst]].obj.show();
			}
			var strand = "";
			if(genes[symbols[css][cst]].strand == "+"){
				strand = ">";
			} else if(genes[symbols[css][cst]].strand == "-"){
				strand = "<";
			}
			cstObj = R.text(x,R_top-sh,genes[symbols[css][cst]].id+strand).attr({font:font_size2_text});
			cstObj[0].style.cursor = "pointer";
			cstObj[0].onclick = function (){
				window.open(ncbi_url+genes[symbols[css][cst]].id);
			};
			cstObj[0].onmouseover = function (){
				//cstObj.attr({font:font_size1_text,"font-weight":"bolder"});
				cstObj[0].style.textDecoration = "underline";
				cstObj.attr({fill:colO_hover2});
			};
			cstObj[0].onmouseout = function (){
				//cstObj.attr({font:font_size2_text,"font-weight":"normal"});
				cstObj[0].style.textDecoration = "";
				cstObj.attr({fill:"#000"});
			};
		}
	}
	function draw_arrow(h1,level,direction){
		var name_length = 100;
		var h = h1 + name_length/2*direction;
		var v = R_top - 20*level;
		var sh = 8;//short height
//		var arrow = R.path("M"+h+","+v+
//					" L"+(h+1*sh*direction)+","+v+
//					" L"+(h+1*sh*direction)+","+(v+sh/2)+
//					" L"+(h+2.5*sh*direction)+","+(v-sh/2)+
//					" L"+(h+1*sh*direction)+","+(v-sh*1.5)+
//					" L"+(h+1*sh*direction)+","+(v-sh)+
//					" L"+h+","+(v-sh)+"Z").attr({fill:colO_hover,"stroke-width":0});					
		if(direction==1){
			var arrow = R.path("M16,1.466C7.973,1.466,1.466,7.973,1.466,16c0,8.027,6.507,14.534,14.534,14.534c8.027,0,14.534-6.507,14.534-14.534C30.534,7.973,24.027,1.466,16,1.466zM13.665,25.725l-3.536-3.539l6.187-6.187l-6.187-6.187l3.536-3.536l9.724,9.723L13.665,25.725z");
		}else {
			var arrow = R.path("M16,30.534c8.027,0,14.534-6.507,14.534-14.534c0-8.027-6.507-14.534-14.534-14.534C7.973,1.466,1.466,7.973,1.466,16C1.466,24.027,7.973,30.534,16,30.534zM18.335,6.276l3.536,3.538l-6.187,6.187l6.187,6.187l-3.536,3.537l-9.723-9.724L18.335,6.276z");
		}
		arrow.attr({fill:colO_hover,stroke:"none"})
		arrow.stop().animate({transform: "t"+ (h-15) +" "+ (v-25) +" s0.7"});
		//			R.path("M"+h+","+v+
		//			" L"+(h+2*sh*direction)+","+(v-sh)+
		//			" L"+h+","+(v-2*sh)+"Z").attr({fill:colO_hover,"stroke-width":0});
	
		(function (obj){
			obj[0].style.cursor = "pointer";
			obj[0].onmouseover = function() {
				obj.attr({fill:colO_hover2});
				R.safari();
			};
			obj[0].onmouseout = function() {
				obj.attr({fill:colO_hover});
				R.safari();
			};
			obj[0].onclick = function() {
				if(level == 0){
					change_transcript(h1,sh,direction);
				}else if(level == 1){
					change_symbol(h1,sh,direction);
				}
			};
		})(arrow);
	}
	function draw_var(natual_pos,name_pos,l,box_width,letter){
		var var_length = 10;
		var var_height = 5;
		var horizontal_offset = R_width-R_right+15+box_width;
		var horizontal_line_offset = horizontal_offset+var_height+var_length;
		var horizontal_text_offset = horizontal_line_offset+50;
		var font_size2_text = "12px \"Trebuchet MS\", Arial, sans-serif";

		var temp = R.set();
		temp.push(R.path("M"+horizontal_offset+","+natual_pos+
			" L"+(horizontal_offset+var_height)+","+(natual_pos+var_height/2)+
			" L"+(horizontal_offset+var_height+var_length)+","+(natual_pos+var_height/2)+
			" L"+(horizontal_offset+var_height+var_length)+","+(natual_pos-var_height/2)+
			" L"+(horizontal_offset+var_height)+","+(natual_pos-var_height/2)+"Z").attr({fill:colO,"stroke-width":0}));

		temp.push(R.path("M"+(horizontal_line_offset)+","+natual_pos
			+" L"+(horizontal_line_offset+10)+","+natual_pos
			+" L"+(horizontal_text_offset-10)+","+name_pos
			+" L"+(horizontal_text_offset)+","+name_pos).attr({stroke:colO,opacity:1}));

		if(letter == "(" || letter == ")"){
			letter = "ASS/DSS loss";
		}else if (letter == "#"){
			letter = "Frame shifting";
		}else if (letter == "^"){
			letter = "Initiator loss";
		}else{
			var texts = letter.split(":");
			if(texts[0].indexOf("$") >= 0 && texts[1].indexOf("$") >= 0){
				letter = "Stop -> Stop";
			}else if(texts[0].indexOf("$") >= 0){
				letter = "Stop loss";
			}else if(texts[1].indexOf("$") >= 0){
				letter = "Stop gain";
			}else {
				letter = texts[0] + "->" + texts[1];
			}
		}

		temp.push(R.text(horizontal_text_offset+2,name_pos,letter).attr({font:font_size2_text,"text-anchor":"start"}));
		return temp;
	}

	function draw_box(from,to,l,box_width){
		var sub_top = map_coord(l,parseInt(from));
		var sub_bot = map_coord(l,parseInt(to)+1);
		if(sub_top == l || sub_bot == 0){
			return null;
		}
		return R.rect(R_width-R_right+10,R_top+sub_top,box_width,sub_bot-sub_top,0).attr({fill:colO,stroke:colO});
	}
	function draw_triangle(from,to,l,strand,box_width){
		var sub_top = map_coord(l,parseInt(from));
		var sub_bot = map_coord(l,parseInt(to)+1);
		if(strand == "+"){
			return R.path("M"+(R_width-R_right+10)+","+(R_top+sub_bot)+" L"+(R_width-R_right+10+box_width)+","+(R_top+sub_bot)+" L"+(R_width-R_right+10+box_width/2)+","+(R_top+sub_bot+box_width)+" Z").attr({fill:colO});
		} else if (strand == "-"){
			return R.path("M"+(R_width-R_right+10)+","+(R_top+sub_top)+" L"+(R_width-R_right+10+box_width)+","+(R_top+sub_top)+" L"+(R_width-R_right+10+box_width/2)+","+(R_top+sub_top-box_width)+" Z").attr({fill:colO});
		} else {
			return null;
		}
	}
	function draw_band(from,to,l,band_width,box_width){
		var sub_top = map_coord(l,parseInt(from));
		var sub_bot = map_coord(l,parseInt(to)+1);
		var slimmer = (box_width-band_width)/2;
		if(sub_top == l || sub_bot == 0){
			return null;
		}
		return R.rect(R_width-R_right+10+slimmer,R_top+sub_top,band_width,sub_bot-sub_top,0).attr({fill:colO,stroke:colO});
	}
	function draw_line(from,to,l,strand,box_width){
		var sub_top = map_coord(l,parseInt(from));
		var sub_bot = map_coord(l,parseInt(to)+1);
		var mid_point_x = R_width-R_right+10+box_width/2;
		if(sub_top == l || sub_bot == 0){
			return null;
		}
		if(strand == "+"){
		//	mid_point_x = R_width-R_right+10+box_width;
			mid_point_x = R_width-R_right+10+box_width/2;
		} else if(strand == "-"){
		//	mid_point_x = R_width-R_right+10;
			mid_point_x = R_width-R_right+10+box_width/2;
		}
		return R.path("M"+(R_width-R_right+10+box_width/2)+","+(R_top+sub_top)+" L"+mid_point_x+","+(R_top+(sub_bot+sub_top)/2)+" L"+(R_width-R_right+10+box_width/2)+","+(R_top+sub_bot));
	}
	function map_coord(l,coord){
		if (coord < current_start){
			return 0;
		} else if (coord > current_end){
			return l;
		} else {
			return (coord-current_start)/(current_end-current_start+1)*l;
		}
	}
}

function if_recombination(recomGB){
	var rcbflag = 0;
	var num = 0;
	var rcbtop = 0;
	var rcbbot = 0;
	var time_t = 5;
	var time_b = 5;
	if(recomGB.length>0){
		for(var i=0; i<recomGB.length; i++){
			if(recomGB[i] != undefined && time_t > 0){
				rcbtop = rcbtop + recomGB[i].value;
				time_t = time_t - 1;
			}
		}
		for(var j=recomGB.length-1; j>0; j--){
			if(recomGB[j] != undefined && time_b > 0){
				rcbbot = rcbbot + recomGB[j].value;
				time_b = time_b - 1;
			}
		}
		if(2*rcbtop == rcbbot){
			rcbflag = 1;//recombination top with blue
		}else if(rcbtop == 2*rcbbot){
			rcbflag = 2;//recombiantion top with green
		}
	}
	return rcbflag;
}

function newif_recombination(recomGB, num){
	var mark = 0;
	var recordpos = [];
	var rep_i;
	var times = 0;
	var rcblocation = [];
	var rcb_i = 0;
	if(recomGB.length>0){
		for(var i=0; i<recomGB.length; i++){
			if(mark == 0){
				mark = recomGB[i].value;
				times = times + 1;
				rep_i = 0;
			}else if(mark == recomGB[i].value){
				times = times + 1;
			}else if(mark != recomGB[i].value){
				recordpos[rep_i] = {};
				recordpos[rep_i].times = times;
				recordpos[rep_i].value = recomGB[i-1].value;
				recordpos[rep_i].pos = recomGB[i-1].pos;
				recordpos[rep_i].nextpos = recomGB[i].pos;
				rep_i = rep_i + 1;
				mark = recomGB[i].value;
				times = 1;
			}
		}
		recordpos[rep_i] = {};
		recordpos[rep_i].times = times;
		recordpos[rep_i].value = recomGB[recomGB.length-1].value;
		recordpos[rep_i].pos = recomGB[recomGB.length-1].pos;
		recordpos[rep_i].nextpos = 0;
		for(var j=0; j<recordpos.length-1; j++){
			if(recordpos[j].times >= num && recordpos[j+1].times >= num && ((recordpos[j].value + recordpos[j+1].value == 7) || (recordpos[j].value + recordpos[j+1].value == 3))){
				rcblocation[rcb_i] = {};
				rcblocation[rcb_i].pos = recordpos[j].pos;
				rcblocation[rcb_i].value = recordpos[j].value;
				rcblocation[rcb_i].nextpos = recordpos[j].nextpos;
				rcb_i = rcb_i +1;
			}else if(j < (recordpos.length-3) && recordpos[j].times >= num && recordpos[j+3].times >= num && ((recordpos[j].value + recordpos[j+3].value == 7) || (recordpos[j].value + recordpos[j+3].value == 3))){
				if(recordpos[j+1].times <= 2 && recordpos[j+2].times <= 2 && (recordpos[j].value == recordpos[j+2].value || recordpos[j+1].value == recordpos[j+3].value)){
					rcblocation[rcb_i] = {};
					rcblocation[rcb_i].pos = recordpos[j].pos;
					rcblocation[rcb_i].value = recordpos[j].value;
					rcblocation[rcb_i].nextpos = recordpos[j].nextpos;
					rcb_i = rcb_i +1;
				}
			}else if(j < (recordpos.length-2) && recordpos[j].times >= num && recordpos[j+2].times >= num && ((recordpos[j].value + recordpos[j+2].value == 7) || (recordpos[j].value + recordpos[j+2].value == 3))){
				if(recordpos[j+1].times <= 2){
					rcblocation[rcb_i] = {};
					rcblocation[rcb_i].pos = recordpos[j].pos;
					rcblocation[rcb_i].value = recordpos[j].value;
					rcblocation[rcb_i].nextpos = recordpos[j].nextpos;
					rcb_i = rcb_i +1;
				}
			}
		}
	}
	return rcblocation;
}

function recomb_loaction(recomGB){
	var record = 0;
	var recom_loc = {};
	recom_loc.pos = -1;
	recom_loc.gap = -1;
	recom_loc.uppos = -1
	for(var i=0; i<recomGB.length; i++){
		if(record == 0){
			record = recomGB[i].value;
		}else if(record != recomGB[i].value){
			recom_loc.pos = recomGB[i].pos
			recom_loc.gap = recomGB[i].pos - recomGB[i-1].pos;
			recom_loc.uppos = recomGB[i-1].pos;
			return recom_loc;
		}
	}
	return recom_loc;
}

function call_siblings_recom(){
	if(reqsibling.readyState == 4) {
		if(reqsibling.status == 200) {
			
			var brwplotCanvas = $(document.getElementById("brwplot"));	
			var topp = brwplotCanvas.position().top + R_top + 195;
			var leftt = brwplotCanvas.position().left + R_left + (R_width-R_left-R_right)/2 + 130;
			$(document.getElementById("sibling_trio_window")).css("top",topp);
			$(document.getElementById("sibling_trio_window")).css("left",leftt);
			$(document.getElementById("sibling_trio_window")).css("display","block");
			
			var RS;
			var RS_height = 425;
			var RS_width = 150;
			var RS_top = 20;
			var RS_left = 25;
			var RS_right = 10;
			var RS_bottom = 20;
			var RS_sremove = null;
			RS = Raphael("sibling_trio", RS_width, RS_height);
	
			var l = RS_height - RS_top - RS_bottom;
			var w = RS_width - RS_left - RS_right;
			var font_size2_text = "10px \"Trebuchet MS\", Arial, sans-serif";
			var axis = RS.path("M"+RS_left+","+RS_top+" L"+RS_left+","+(RS_height-RS_bottom)).attr({stroke:colO,"stroke-width":1});
			var cali = [];
			/*var warningtext1 = RS.text(0,RS_top-40,"The requested scope is too large.").attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colD});
			var warningtext2 = RS.text(0,RS_top-25,"Please wait ...").attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colD});
			var warning_sibling_set = RS.set();
			warning_sibling_set.push(warningtext1, warningtext2);
			if(current_end - current_start < 3000000){
				warning_sibling_set.hide();
			}else{
				warning_sibling_set.show();
			}*/
			var t_bot, t_top;
			t_top = RS.text(0,RS_top-10,current_chr+":"+current_start).attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colO});
			t_bot = RS.text(0,RS_height-RS_bottom+10,current_chr+":"+current_end).attr({font: font_size2_text ,opacity:1,"text-anchor":"start"}).attr({fill: colO});
			//axis_set.push(axis, t_bot, t_top);
			for(var ca = 0 ; ca <= 50 ; ca++){
				if(ca%10 == 0){
					cali[ca] = RS.path("M"+(RS_left-20)+","+(RS_top+ca/50*l)+" L"+RS_left+","+(RS_top+ca/50*l)).attr({stroke:colO,"stroke-width":1});
				}else{
					cali[ca] = RS.path("M"+(RS_left-10)+","+(RS_top+ca/50*l)+" L"+RS_left+","+(RS_top+ca/50*l)).attr({stroke:colO,"stroke-width":1});
				}
				//axis_set.push(cali[ca]);
			}
			var bin_size = 5;
			var pixel_per_variant = 1;
			var l = RS_height - RS_top - RS_bottom;
			var w = RS_width - RS_left - RS_right;
			RS.rect(RS_left+w/2-15,RS_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
			RS.rect(RS_left+w/2+5,RS_top,10,l,5).attr({fill:colO_hover,stroke:colO_hover});
			
			var redbins = [];
			redbins[0] = [];
			redbins[1] = [];
			var greenbins = [];
			greenbins[0] = [];
			greenbins[1] = [];
			var bluebins = [];
			bluebins[0] = [];
			bluebins[1] = [];

			var sib_bin_set = RS.set();
			var sib_recom_GB = [];
			sib_recom_GB[0] = [];
			sib_recom_GB[1] = [];


			if(!outScopeflag){
				var vsNodes = reqsibling.responseXML.getElementsByTagName(xmlTagVariants);
				for(var i = 0 ; i < vsNodes.length ; i++){
					var vNodes = vsNodes[i].getElementsByTagName(xmlTagVariant);
					var vs_id = vsNodes[i].getAttribute(xmlAttributeId);
					var vPointer = 0;
					for(var j = 0 ; j < vNodes.length ; j++){
						var from = vNodes[j].getElementsByTagName(xmlTagFrom)[0].childNodes[0].nodeValue;
						var to = vNodes[j].getElementsByTagName(xmlTagTo)[0].childNodes[0].nodeValue;
						var id = vNodes[j].getAttribute(xmlAttributeId);
						var type = vNodes[j].getAttribute(xmlAttributeType);

						if(id == "."){
							id = "Unnamed "+type;
						}
						while(variants[vPointer].from < from){
							vPointer++;
						}
						if(from <= variants[vPointer].from && to >= variants[vPointer].to && id == variants[vPointer].id){
							if(vs_id == individuals[csi].fid){
								variants[vPointer].paternal = "Y";
								variants[vPointer].maternal = "N";
								change_variant_color(vPointer,colP);
								var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
								if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
									if(bluebins[0][Math.floor(natual_pos)] == undefined){
										bluebins[0][Math.floor(natual_pos)] = 1;
									}else{
										bluebins[0][Math.floor(natual_pos)] ++;
									}
								}
								if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
									if(bluebins[1][Math.floor(natual_pos)] == undefined){
										bluebins[1][Math.floor(natual_pos)] = 1;
									}else{
										bluebins[1][Math.floor(natual_pos)] ++;
									}
								}
							}else if(vs_id == individuals[csi].mid){
								variants[vPointer].paternal = "N"
									variants[vPointer].maternal = "Y";
								change_variant_color(vPointer,colM);
								var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*(l/bin_size);
								if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
									if(greenbins[0][Math.floor(natual_pos)] == undefined){
										greenbins[0][Math.floor(natual_pos)] = 1;
									}else{
										greenbins[0][Math.floor(natual_pos)] ++;
									}
								}
								if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
									if(greenbins[1][Math.floor(natual_pos)] == undefined){
										greenbins[1][Math.floor(natual_pos)] = 1;
									}else{
										greenbins[1][Math.floor(natual_pos)] ++;
									}
								}
							}else if(vs_id == csi){
								variants[vPointer].paternal = "N";
								variants[vPointer].maternal = "N";
								change_variant_color(vPointer,colD);

								if(variants.length > l/3){
									var natual_pos = ((variants[vPointer].to+variants[vPointer].from)/2 - current_start)/(current_end - current_start + 1)*l;
									if(variants[vPointer].genotypes[0] != undefined && variants[vPointer].genotypes[0] != "0"){
										redbins[0][j] = natual_pos;
									}
									if(variants[vPointer].genotypes[1] != undefined && variants[vPointer].genotypes[1] != "0"){
										redbins[1][j] = natual_pos;
									}
								}
							}
						}
					}
				}
			}else{
				var vlNodes = reqsibling.responseXML.getElementsByTagName(xmlTagValues);
				for(var i=0; i<vlNodes.length; i++){
					switch(vlNodes[i].getAttribute(xmlAttributeId)){
						case (families[csi].fid + "L"):
							bluebins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
						case (families[csi].fid + "R"):
							bluebins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
						case (families[csi].mid + "L"):
							greenbins[0] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
						case (families[csi].mid + "R"):
							greenbins[1] = vlNodes[i].getElementsByTagName(xmlTagValueList)[0].firstChild.nodeValue.split(",");
							break;
					}
				}
			}
			if(outScopeflag || variants.length >= l/3){
				var ratio = 0;
				var color = [];
				color[0] = colP;
				color[1] = "#83A3CD";
				color[2] = "#C1D1E6";
				color[3] = "#FFFFFF";
				color[4] = "#B3DCC4";
				color[5] = "#66BA8A";
				color[6] = colM;
				color[7] = colD;
				var colBorder = colO_hover2;
				var sw = 1;
				var re_i = 0;
				var if_firstB = [];
				var if_firstG = [];
				for(var i=0; i<bins[0].length; i++){
					if(greenbins[0][i] == undefined){
						greenbins[0][i] = 0;
					}
					if(bluebins[0][i] == undefined){
						bluebins[0][i] = 0;
					}
					if(Math.abs(bluebins[0][i]) > 0 || Math.abs(greenbins[0][i]) > 0){
						if(Math.abs(bluebins[0][i]) >= Math.abs(greenbins[0][i])){
							if(if_firstB[0] == undefined){
								if(bluebins[0][i] > 0){
									if_firstB[0] = 1;
								}else{
									if_firstB[0] = -1;
								}
							}	
							if(bluebins[0][i] > 0){
								if(if_firstB[0] > 0){
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 1;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 2;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstB[0] > 0){
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 2;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 1;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}
							}
						}else if(Math.abs(bluebins[0][i]) < Math.abs(greenbins[0][i])){
							if(if_firstG[0] == undefined){
								if(greenbins[0][i] > 0){
									if_firstG[0] = 1;
								}else{
									if_firstG[0] = -1;
								}
							}
							if(greenbins[0][i] > 0){
								if(if_firstG[0] > 0){
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 4;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 3;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstG[0] > 0){
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 3;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[0][re_i]={};
									sib_recom_GB[0][re_i].value = 4;
									sib_recom_GB[0][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2-40+6,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}
							}
						}
					}
				}
				

				re_i = 0;
				for(var i=0; i<bins[1].length; i++){
					if(greenbins[1][i] == undefined){
						greenbins[1][i] = 0;
					}
					if(bluebins[1][i] == undefined){
						bluebins[1][i] = 0;
					}

					if(Math.abs(bluebins[1][i]) > 0 || Math.abs(greenbins[1][i]) > 0){
						if(Math.abs(bluebins[1][i]) >= Math.abs(greenbins[1][i])){
							if(if_firstB[0] == undefined){
								if(bluebins[1][i] > 0){
									if_firstB[0] = 1;
								}else{
									if_firstB[0] = -1;
								}
							}	
							if(bluebins[1][i] > 0){
								if(if_firstB[0] > 0){
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 1;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 2;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstB[0] > 0){
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 2;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[1],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 1;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[0],stroke:colBorder,"stroke-width":1}));
								}
							}
						}else if(Math.abs(bluebins[1][i]) < Math.abs(greenbins[1][i])){
							if(if_firstG[0] == undefined){
								if(greenbins[1][i] > 0){
									if_firstG[0] = 1;
								}else{
									if_firstG[0] = -1;
								}
							}
							if(greenbins[1][i] > 0){
								if(if_firstG[0] > 0){
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 4;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 3;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}
							}else{
								if(if_firstG[0] > 0){
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 3;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[5],stroke:colBorder,"stroke-width":1}));
								}else{
									sib_recom_GB[1][re_i]={};
									sib_recom_GB[1][re_i].value = 4;
									sib_recom_GB[1][re_i].pos = i;
									re_i = re_i+1;
									sib_bin_set.push(RS.rect(RS_left+w/2+40-16,i*bin_size+RS_top,10,bin_size,0).attr({fill:color[6],stroke:colBorder,"stroke-width":1}));
								}
							}
						}
					}
				}
				for(var j=0; j<redbins[0].length; j++){
					if(redbins[0][j] != undefined){
						sib_bin_set.push(RS.path("M "+ (RS_left+w/2-40+6) + "," + (redbins[0][j]+RS_top) + " L" + (RS_left+w/2-40+16) + "," + (redbins[0][j]+RS_top)).attr({stroke:colD,"stroke-width":1}).toFront());
					}
				}
				for(var j=0; j<redbins[1].length; j++){
					if(redbins[1][j] != undefined){
						sib_bin_set.push(RS.path("M "+ (RS_left+w/2+40-16) + "," + (redbins[1][j]+RS_top) + " L" + (RS_left+w/2+40-6) + "," + (redbins[1][j]+RS_top)).attr({stroke:colD,"stroke-width":1}).toFront());
					}
				}
				var chr_left;
				var chr_right;
				if(current_end - current_start <= 20000000){
					chr_left = newif_recombination(sib_recom_GB[0],5);
					chr_right = newif_recombination(sib_recom_GB[1],5);
				}else{
					chr_left = newif_recombination(sib_recom_GB[0],3);
					chr_right = newif_recombination(sib_recom_GB[1],3);
				}
				var rcbicon1 = RS.path("M26.711,14.086L16.914,4.29c-0.778-0.778-2.051-0.778-2.829,0L4.29,14.086c-0.778,0.778-0.778,2.05,0,2.829l9.796,9.796c0.778,0.777,2.051,0.777,2.829,0l9.797-9.797C27.488,16.136,27.488,14.864,26.711,14.086z").attr({fill: "#ccc", stroke: "none"});
				var rcbicon2 = RS.path("M14.702,8.981c0.22-0.238,0.501-0.357,0.844-0.357s0.624,0.118,0.844,0.353c0.221,0.235,0.33,0.531,0.33,0.885c0,0.306-0.101,1.333-0.303,3.082c-0.201,1.749-0.379,3.439-0.531,5.072H15.17c-0.135-1.633-0.301-3.323-0.5-5.072c-0.198-1.749-0.298-2.776-0.298-3.082C14.372,9.513,14.482,9.22,14.702,8.981z").attr({fill: "#000", stroke: "none"});
				var rcbicon3 = RS.path("M16.431,21.799c-0.247,0.241-0.542,0.362-0.885,0.362s-0.638-0.121-0.885-0.362c-0.248-0.241-0.372-0.533-0.372-0.876s0.124-0.638,0.372-0.885c0.247-0.248,0.542-0.372,0.885-0.372s0.638,0.124,0.885,0.372c0.248,0.247,0.372,0.542,0.372,0.885S16.679,21.558,16.431,21.799z").attr({fill: "#000", stroke: "none"});
				rcbicon1.hide();
				rcbicon2.hide();
				rcbicon3.hide();
				var sib_recom_set = RS.set();
				var icon_set = RS.set();
				sib_recom_set.push(rcbicon1,rcbicon2,rcbicon3);
				if(chr_left.length != 0){
					for(var i=0; i<chr_left.length; i++){
						var clone1 = rcbicon1.clone();
						var clone2 = rcbicon2.clone();
						var clone3 = rcbicon3.clone();
						clone1.animate({transform:"t"+(RS_left+w/2-40+16)+" "+(RS_top+chr_left[i].pos*bin_size-11)});
						clone2.animate({transform:"t"+(RS_left+w/2-40+16)+" "+(RS_top+chr_left[i].pos*bin_size-11)});
						clone3.animate({transform:"t"+(RS_left+w/2-40+16)+" "+(RS_top+chr_left[i].pos*bin_size-11)});
						sib_recom_set.push(clone1,clone2,clone3);
						icon_set.push(clone2,clone3);
					}
				}
				if(chr_right.length != 0){
					for(var i=0; i<chr_right.length; i++){
						var clone1 = rcbicon1.clone();
						var clone2 = rcbicon2.clone();
						var clone3 = rcbicon3.clone();
						clone1.animate({transform:"t"+(RS_left+w/2+40-46)+" "+(RS_top+chr_right[i].pos*bin_size-11)});
						clone2.animate({transform:"t"+(RS_left+w/2+40-46)+" "+(RS_top+chr_right[i].pos*bin_size-11)});
						clone3.animate({transform:"t"+(RS_left+w/2+40-46)+" "+(RS_top+chr_right[i].pos*bin_size-11)});
						sib_recom_set.push(clone1,clone2,clone3);
						icon_set.push(clone2,clone3);
					}
				}

				var anim_0 = Raphael.animation({
					100:{
						fill:"#000"
					},
					50:{
						fill:"#ccc"
					}},1000);
				icon_set.animate(anim_0.repeat("Infinity"));
			}
			document.body.addEventListener("mousedown", mousedownOutsideRepeatTooltip, false);
			function mousedownOutsideRepeatTooltip(evt){
				evt = evt || window.event;
				var eventTarget = evt.target || evt.srcElement;
				var flag=0;
				while(eventTarget){
					if(eventTarget==document.getElementById("sibling_trio_window")){
						flag=1;
						break;
					}else{
						eventTarget = eventTarget.parentNode;
					}
				}
				if(flag==0) {
					document.getElementById("sibling_trio_window").style.display = "none";
				}
			}
		}
	}
}
function close_siblings_recom(){
	document.getElementById("sibling_trio_window").style.display = "none";
}
