#include "Cello.h"

// author: Sukriti Bhattacharya
// canopy clustering algorithm using cello
var read_from_file(var args)
{
	var ct = get(args, $I(0)); // contig table
	var data = new(List, Float);
	var null = new(List, Float);
	var length = new(Array, Float); // List contaning the length of the contigs, to be sorted
	var srt_ct = new(Tuple);  // < ct, length>	
	FILE *fp; // file pointer

	var contig = new(String, $S("contig_"));
	var f0 =new(Int, $I(0)); // contig id after "_"

	var f1 =new(Float);
	var f2 =new(Float);var f3 =new(Float);var f4 =new(Float);var f5 =new(Float);
	var f6 =new(Float);var f7 =new(Float);var f8 =new(Float);var f9 =new(Float);

	fp = fopen("cstr_wh.csv","r"); // open the file

	char buff[255]; // each line of the file, i.e. each contig
	char itoa[20];  // convert f0 to string to form "contig_f0" 
	var str = copy(contig); // str == "contig_"
	while(!feof(fp)){
		 fgets(buff, 200, fp);
		 var line = new(String, $S(buff));
		 scan_from(line, 0, "contig_%$,%$,%$,%$,%$,%$,%$,%$,%$,%$", f0,f1,f2,f3,f4,f5,f6,f7,f8,f9);  // format string
		 push(data, f1); push(length, f1);
		 push(data, f2);push(data, f3);push(data, f4);push(data, f5);push(data, f6);
		 push(data, f7);push(data, f8);push(data, f9);
		 int value = c_int(f0);
		 snprintf(itoa, 10, "%d", value);
		 var  id = new(String, $S(itoa));
		 concat(str,id); // "contig_" || "f0"
		 set(ct, str, data);
		 str = copy(contig);
		 data = copy(null);
	}
	fclose(fp);  // close the file
	sort_by(length,gt); // sort by descending order of length
	push(srt_ct, ct);
	push(srt_ct, length);
	return srt_ct;
}

var sort_contig(var args) // associate le
{
	var ct = get(args, $I(0));
	var len = get(args, $I(1));
	var contig_id = new(List, String);
	var ct_id = new(Tuple);
	foreach (l in len){
		foreach (i in ct){	
			if (mem(get(ct, i), l)){
				push(contig_id, i);
				rem(get(ct, i), l);
			}	
		}
	}
	push(ct_id, ct);
	push(ct_id, contig_id);
	return ct_id;	
}

var calculate_pcc(var args)
{
    var x = get(args,$I(0));
    var y = get(args,$I(1));
    var sum_xy = new(Float, $F(0.0));
    var sum_x = new(Float, $F(0.0));
    var sum_y = new(Float, $F(0.0));
    var sum_x2 = new(Float, $F(0.0));
    var sum_y2 = new(Float, $F(0.0));
    var n = copy($F(len(x)));
    var r = new(Float);	
    foreach (pair in zip(x,y)){
	    var xi = get(pair, $I(0));
	    var yi = get(pair, $I(1));
	    var xy = $F(c_float(xi) * c_float(yi));
	    var x2 = $F(pow(c_float(xi) ,2.0));
	    var y2 = $F(pow(c_float(yi) ,2.0));
            assign(sum_x, $F(c_float(sum_x) + c_float(xi)));
	    assign(sum_y, $F(c_float(sum_y) + c_float(yi)));
	    assign(sum_x2, $F(c_float(sum_x2) + c_float(x2)));
            assign(sum_y2, $F(c_float(sum_y2) + c_float(y2)));
            assign(sum_xy, $F(c_float(sum_xy) + c_float(xy)));
    }
    var a = $F(c_float(n) * c_float(sum_xy) - c_float(sum_x) * c_float(sum_y));
    var b1 = $F((c_float(n) * c_float(sum_x2)- pow(c_float(sum_x), 2.0)));
    var b2 = $F((c_float(n) * c_float(sum_y2)- pow(c_float(sum_y), 2.0)));
    var b = $F(sqrt(c_float(b1)*c_float(b2)));
    assign(r, $F(c_float(a)/c_float(b)));
    return r;
}


var cluster_contigs(var args)
{
 var ct = get(args, $I(0));
 var id = get(args, $I(1));
 var k = new(Int, $I(0));
 var cluster = new(Table, String, List);
 var t = new(Float, $F(0.90));
 foreach(i in id){
	assign(k, $I(c_int(k)+1));
	var elements = new(List, String);
	foreach(j in slice(id, k, $I(len(id)))){
		var r = call($(Function, calculate_pcc), get(ct,i), get(ct,j));
		if (gt(r, t)){
			push(elements, j);
			rem(id, j);
		}
	}
	set(cluster, i, elements);
 }	
 
 return cluster;	
}

var merge(var args)
{
 	var cluster_old = get(args, $I(0));
	var contig_list = get(args, $I(1));
	var l = new(List, String);
	foreach(i in contig_list){
		var list = get(cluster_old, i);
		concat(l, list);
	}
	return l;
		
}

var final_cluster(var args)
{
	var cluster_centroid = get(args, $I(1));
	var cluster_old = get(args, $I(2));
       	var cluster_new = new(Table, String, List);
	
	foreach(i in cluster_centroid){
			var n = copy($I(len(get(cluster_centroid,i))));	
			if(ge(n, $I(2))){
				var ml = call($(Function, merge), cluster_old, get(cluster_centroid,i));	
				set(cluster_new, i, ml);
				continue;
			}
			set(cluster_new, i, get(cluster_old, i));
	}
		
	
        return cluster_new;
}


var cluster_centroids(var args)  
{
	var centroids = get(args, $I(0));   
	var cluster_old = get(args, $I(1)); 
	
	var id = new(List, String);
	foreach(i in centroids)
	{
		push(id,i);
	}	
	var cluster_centroid = call($(Function, cluster_contigs), centroids, id);
	var cluster_final = call($(Function, final_cluster), centroids, cluster_centroid, cluster_old);
	return cluster_final;
}



var add(var args)
{
 var l = new(List, Float);
 foreach(i in zip (get(args, $I(0)), get(args, $I(1))))
  {
   var sum = new(Float);
   assign(sum, $F(c_float(get(i,$I(0))) + c_float(get(i,$I(1)))));
   push(l,sum);
  }
 return l;
}

var mean(var args)
{
   var l =new(List, Float);
   var n = get(args,$I(1));
   foreach (i in get(args, $I(0))){
     var m = new(Float);
     assign(m,$F(c_float(i)/c_float(n))); //devide each eliment of list by the total number of points in the cluster 
     push(l,m);     
   }
     return l;
}


var calculate_mean(var args)
{
  var contig = get(args, $I(0));
  var centroid = get(args, $I(1));
  var cluster = get(args, $I(2));
  var n = new(Float);
  foreach (i in cluster){
	var l1 = new(List, Float); 	
	var l0 = new(List, Float, $F(0.0), $F(0.0), $F(0.0),$F(0.0), $F(0.0), $F(0.0), $F(0.0), $F(0.0));
	var n = copy($F(len(get(cluster,i))));  // number of contigs in each cluster
	if (ge(n, $F(2.0))){	
		foreach (j in contig){     
			if (mem(get(cluster, i), j)) // if the contig belogs to the cluster
				{
					l1 = copy(call($(Function, add), l0, get(contig,j))); // pairwise addition
					l0 = copy(l1);
        		        } 
		}	
   	
 		var m = call($(Function, mean), l1, n);
  		set(centroid, i, m);
	}
 }
return centroid;
}


var claculate_centroid(var args)
{
	var ct = get(args, $I(0));
	var clusters = get(args, $I(1));	
	var l = new(Int);
	var centroids = new(Table, String, List);    
	centroids = call($(Function, calculate_mean), ct, centroids, clusters);
	return centroids;
}

int main(int argc, char** argv) {
	var contigs = new(Table, String, List);     // contigs <id, pointer to the length, coordinates>
	var contig = call($(Function, read_from_file), contigs); //tuple < contig_table, sorted length array>
	var contig_table = get(contig, $I(0)); // contig table with length field 
	var l = get(contig, $I(1)); // len: sorted length array
	
	var ct_id = call($(Function, sort_contig), contig_table, l); // sorted contig list, based on descending order of length
	
	var ct = get(ct_id, $I(0)); // contig list without length field
	var id = get(ct_id, $I(1)); // contig id, sorted by their length
	var cluster = call($(Function, cluster_contigs), ct, id); // claculate the first cluster based on pcc on sorted contigs
	
	var centroids = call($(Function, claculate_centroid), ct, cluster); // claculate centroid by claculating mean on each cluster 		elements
	var cluster_f = call($(Function, cluster_centroids), centroids, cluster); // cluster the centroids based on PCC, yields, the final cluster

show(cluster_f);
return 1;
}


