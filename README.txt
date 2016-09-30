
 ++++++++ MEMSA v 1.0 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 This program can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software 
 Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed by the author in the 
 hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE. If you use this program please cite 
 
 [1] J. Iacovacci, Z. Wu and G. Bianconi, Mesoscopic structures reveal the network between the layers of multiplex datasets,
 Physical Review E 92 (4), 2015  

 (c) Jacopo Iacovacci (mriacovacci@hotmail.it) 
     Ginestra Bianconi (ginestra.bianconi@gmail.com)
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MEMSA (MEsoscopic Multiplex Structure Analysis) is a program developed to extract information from the mesoscopic structures of multiplex
networks induced by node features. It can be used to asses the similarity of layers of multiplex networks with repect to a generic node features 
depending on the network architecture (for example the community assignement). This analysis can be used to extract a network between the layers 
of multiplex network datasets (see [1]).  

The program takes in input the following sequence 

 [N] [edge_list_layer1.txt] [feature_layer1.txt] [edge_list_layer2.txt] [feature_layer2.txt] [k-classes] [stirling] [Nr] 

  N                    -- integer value corresponding to the number of nodes
 
  edge_list_layer1.txt -- file containing the egde list respect to layer 1 
                          (nodes ID labelled from 0 to N-1) 


  feature_layer1.txt   -- file containing the integer values of the features of corresponding 
                          the nodes in layer 1 
                          (nodes ID labelled from 0 to N-1) 


  edge_list_layer2.txt -- file containing the egde list respect to layer 2 
                          (nodes ID labelled from 0 to N-1) 


  feature_layer2.txt   -- file containing the integer values of the features of corresponding 
                          the nodes in layer 2 
                          (nodes ID labelled from 0 to N-1) 


  k-classes            -- set classes of nodes respect to their degree 
                          (0 = off, 1 = on).
 

  stirling             -- set the Stirling approximation to compute entropies 
                          (0 = off, 1 = on)
 

  Nr                   -- numbers of reshuffling to compute the Z-score function ThetaS
                          


and it provides the following output structure 

Theta(1,1)  -- element (1,1) of the matrix Theta 
               (entropy Z-score indicating the significance of the node feature in 
               layer 1 with respect to the network structure of layer 1)


  
Theta(1,2)  -- element (1,2) of the matrix Theta 
               (entropy Z-score indicating the significance of the node feature in 
               layer 2 with respect to the network structure of layer 1)
 

Theta(2,2)  -- element (2,2) of the matrix Theta 
               (entropy Z-score indicating the significance of the node feature in 
               layer 2 with respect to the network structure of layer 2)


Theta(2,1)  -- element (2,1) of the matrix Theta 
               (entropy Z-score indicating the significance of the node feature in 
               layer 1 with respect to the network structure of layer 2)



ThetaS score -- indicates the similarity between layer 1 and layer 2 with
                respect to the node feature 


it also provides informations regarding the dimensionality of the block model 
used to evaluate the etropies    

dimensionality of layer 1 -- numebr of blocks in layer 1  
  
dimensionality of layer 2 -- number of blocks in layer 2


For a more detailed description of the measures see

 [1] J. Iacovacci, Z. Wu and G. Bianconi, Mesoscopic structures reveal the network between the layers of multiplex datasets,
 Physical Review E 92 (4), 2015  

********************************************************************************************
QUICK START
********************************************************************************************


Inside the folder MEMSA you find 


1) thetaS_analysis.c   (source code)
2) Small_Network       (a folder containing a two-layer multiplex network dataset of 61 nodes
                        where the node feature corresponds to the community index of the nodes) 


To test the code 

1) open the MEMSA folder 
2) extract the compressed folder Small_Network.zip 
3) open the terminal inside the MEMSA folder and type 
   the following instructions


>> gcc -g thetaS_analysis.c -o thetaS_analysis -lm

>> cp thetaS_analysis Small_Network/thetaS_analysis

>> ./thetaS_analysis 61 edge_list_1.txt comm_list_1.txt edge_list_2.txt comm_list_2.txt 0 0 500


and you will get an output of the form 


dimensionality of layer 1 = 135  
dimensionality of layer 2 = 60

Theta(1,1) = -3.63673 
Theta(1,2) = -2.87483 
Theta(2,1) = -2.34605 
Theta(2,2) = -2.80258 

ThetaS score = 0.813803 



