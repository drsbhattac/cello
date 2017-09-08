#include "Cello.h"
// author: Sukriti Bhattacharya
// kmean clustering algorithm using cello

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
  var new_centroid = new(Table, String,List);
  new_centroid = copy(centroid);
 
  foreach (i in cluster){
	var l1 = new(List, Float); 	
	var l0 = new(List, Float, $F(0.0), $F(0.0), $F(0.0),$F(0.0), $F(0.0), $F(0.0));
  	foreach (j in contig){     
    		if (mem(get(cluster, i), j)) // if the contig belogs to the cluster
		{
			l1 = copy(call($(Function, add), l0, get(contig,j))); // pairwise addition
			l0 = copy(l1);
                } 
		
   	}
  var n = copy($F(len(get(cluster,i))));  // number of contigs in each cluster
  if (neq(n,$F(0.0))){			  // only for non empty clusters
  	var m = call($(Function, mean), l1, n);
	set(new_centroid, i, m);
  }  
}
 return new_centroid;
}

var remove_element(var args)
{
  var clusters = get(args,$I(0)); 
  foreach(i in clusters){
	if(mem(get(clusters, i), get(args,$I(1))))
	   rem(get(clusters, i), get(args,$I(1)));
  }
  	
 return clusters;
}

var add_to_cluster(var args)
{
  var dist = get(args, $I(0));
  var clusters = copy(get(args,$I(1)));
  foreach(i in dist)
  {
    var key = new(String);
    key = get(i, $I(1));
    if (!(mem(get(clusters,key), get(i,$I(0))))) // if the contif not belongs to List then only add it	
	{	
		// then the contig must be assigned to some cluster. remove it first, before adding.   		
		clusters = call($(Function, remove_element), clusters, get(i,$I(0)));  
		push(get(clusters,key), get(i,$I(0))); // add contig to the cluster
	}
  }
   return clusters;
}


var find_min(var args)
{
 var t= new(Tuple);
 assign(t,get(args,$I(0)));
 var mt = new(Tuple); 
 var min = new(Float);
 assign(min, get(get(t,$I(0)),$I(2)));
 assign(mt,get(t,$I(0)));
 foreach (s in slice(t, $I(1),$I(len(t)))) 
 {
    if(lt(get(s,$I(2)),min)){ 
 	assign(min,get(s,$I(2)));
        assign(mt,s);
    }
}
return mt;
}


var calculate_distance(var args)
{
 var x = get(args, $I(0));
 var y = get(args, $I(1));
 var sum = $(Float, 0.0);
 
 foreach (pair in zip(x, y)) {
       var k = new(Float);
       assign(k, $F(pow((c_float(get(pair,$I(0)))-c_float(get(pair,$I(1)))),2.0)));
       assign(sum, $F(c_float(k)+c_float(sum)));
 }
   var d = new(Float);  
   assign(d,$F(sqrt(c_float(sum))));
   return d;
}


var distance (var args)
{
  var contigs = get(args,$I(0));
  var centroids = get(args,$I(1));
  var min_dist = new(Tuple);
  foreach (i in contigs){ 
	var dist = new(Tuple);
	var mt = new(Tuple);
	foreach(j in centroids){  
		//print("%$", get(centroids, j));	
        	//print("\n");
   		var d = call($(Function, calculate_distance), get(contigs, i), get(centroids, j));
        	// <i, j, distant>
		var z = new(Tuple);
		push(z,i);
        	push(z,j); 
        	push(z,d);
		
		//<<i, j1, distance>, ...<i, j_i, distance>
		push(dist,z);
        	  	 }	
	mt = call($(Function, find_min), dist); // calculate the minimum distance 
    	push(min_dist,mt);
  }
  return min_dist; 
}

var compare(var args)
{
 var old = get(args,$I(0));
 var new = get(args,$I(1));
 var yes = new(Int);
 var no = new(Int);
 assign(yes,$I(1));
 assign(no,$I(0));
 foreach(pair in zip(old,new)){
  var v1 = get(old,get(pair, $I(0))); 
  var v2 = get(new,get(pair, $I(0)));  
  if (neq(v1,v2))
    return no;
 }
 return yes; 
}




int main(int argc, char** argv) {


var contigs = new(Table, String, List);     // contigs <id, pointer to the coordinates>
var centroids = new(Table, String, List);   // centroids <id, pointer to the co ordinates>
var clusters = new(Table, String, List);   // clusters <centroid, list of contigs>

// config information
var x1 = new(List, Float, $F(3.1), $F(5.1), $F(7.1),$F(9.1), $F(11.1), $F(13.0)); 
var x2 = new(List, Float, $F(3.1), $F(5.1), $F(7.1),$F(9.1), $F(10.1), $F(12.1));


// centroid information
var c1 = new(List, Float, $F(2.0), $F(4.0), $F(6.0),$F(8.0), $F(10.0), $F(12.0)); 
var c2 = new(List, Float, $F(3.1), $F(5.1), $F(7.1),$F(9.1), $F(11.1), $F(13.1));

// first cluster elements
var cl1 = new(List, String);
var cl2 = new(List, String);


// populate the data table
set(contigs, $S("contig_1"),   x1); 
set(contigs, $S("contig_2"),   x2);
//set(contigs, $S("contig_3"),   c1); 
//set(contigs, $S("contig_4"),   c2);


// populate the centroid table
set(centroids, $S("contig_3"),   c1); 
set(centroids, $S("contig_4"),   c2);


// add the first eliments of the each cluster
set(clusters, $S("contig_3"), cl1);
set(clusters, $S("contig_4"), cl2); 
var match = new(Int);
do{
	var dist = call($(Function, distance), contigs, centroids); // distance calls 2 functions, calculate_distance and find_min
	var new_clusters = call($(Function, add_to_cluster), dist, clusters); // assign points to clusters
  	var new_centroids = new(Table, String, List); 
	new_centroids = call($(Function, calculate_mean), contigs, centroids, new_clusters);
	match = call($(Function, compare), centroids, new_centroids);
	show(match);	
	centroids = copy(new_centroids);
	clusters = copy(new_clusters);	
}while(neq(match,$I(1)));
show(clusters);
return 1;
}
