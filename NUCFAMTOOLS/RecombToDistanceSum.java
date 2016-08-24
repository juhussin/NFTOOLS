//
//  RecombToDistanceSum.java
//  LVTOOLS
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Summarize recombination events from Recomb to summary files (.females, .males, .parents) taking into account distances from centromere/telomere
//

//package LVTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class RecombToDistanceSum {
	
	
	public static void main(String[] args)
	{
		
		//
		//ARGUMENT TREATMENT
		//
		
		int minArg = 2;
		
		if (args.length < minArg) exit(" usage : java RecombToDistanceSum chrlist.files fams.txt [-mb -rel x -abs x max t/c -o outfileprefix -stop -group file -s bef/aft]");
		

		//Arguments
		
		String filelist = args[0];	//results from recomb
		String famfile = args[1];	//families : famId \t father \t mother \t child1 child2 child3 ... \t ages
		String outfile = filelist;	//outfile
		boolean mb = false;
		
		//default
		boolean relative = true;
		double x = 0.25;
		int abs_max = -1;
		char abs_ref = '?';
		boolean stop = false;
		
		boolean already = false;
		boolean sign = false;
		String signtag = "";
		HashSet Group = null;
		
		if(args.length > minArg)
		{
			int index = minArg;
			while(index < args.length)
			{
				if(args[index].equals("-o"))
				{
					outfile = args[index+1];
					index += 2;
				}
				else if(args[index].equals("-mb"))
				{
					mb = true;
					index += 1;
				}
				else if(args[index].equals("-stop"))
				{
					stop = true;
					index += 1;
				}
				else if(args[index].equals("-s"))
				{
					sign = true;
					signtag = args[index+1];
					index += 2;
				}
				else if(args[index].equals("-group"))
				{
					Group = new HashSet();
					try
					{
						BufferedReader BR = new BufferedReader( new FileReader(args[index+1]));
						String l = BR.readLine(); // assuming that there is no legend!!!!!
						while(l!=null)
						{
							String[] lT = l.split("\t");
							Group.add(lT[0]);
							l = BR.readLine();
						}
						BR.close();
					}
					catch(Exception e)
					{
						e.printStackTrace();
						exit("File reading problem");
					}
					
					
					index += 2;
				}

				else if(args[index].equals("-rel"))
				{
					if (already) exit("Specify only one of the 2 options : [-rel x] or [-abs x max t/c]");
					relative = true;
					x = Double.parseDouble(args[index+1]);
					Double a = new Double(1/x);
					if(1/x != a.intValue()) 
					{
						exit("For -rel : specify x such as 1/x is an integer : "+ (1/x)+" "+a.intValue());
					}
					
					index += 2;
					already = true;
				}
				else if(args[index].equals("-abs"))
				{
					exit("-abs option is unavailable for the moment");
					/*
					if (already) exit("Specify only one of the 2 options : [-rel x] or [-abs x max t/c]");
					relative = false;
					x = Double.parseDouble(args[index+1]);
					abs_max = -1;
					abs_ref = args[index+3].charAt(0);
					if(abs_ref!='t' && abs_ref!='c')
					{
						exit("For -abs : specify t for telomere or c for chromosome : -abs x max t/c");
					}
					if(abs_max%x != 0) 
					{
						exit("For -abs : specify x such as max/x is an integer");
					}
					index += 4;
					already = true;
					*/
				}
				else
				{
					exit("Invalid option : "+ args[index]);
				
				}
			}
			
		
		}
		
		if(relative)
		{
			System.err.println("\nRelative distance will be computed with x = "+x);
		}
		else	//does not work yet!!
		{
			System.err.print("\nAbsolute distance will be computed with x = "+x+" until a maximum distance of "+ abs_max +" from ");
			if (abs_ref == 't') System.err.println("telomeres.");
			if (abs_ref == 'c') System.err.println("centromere.");
		}
		
		//
		//END OF ARGUMENT TREATMENT
		//
		
		
		
		//
		// GETTING THE DATA
		//
		
		//1- Number of files (bins of distances) + binsname table
		
		int nbbins = -1;
		String[] binsname = null;
		
		if(relative)
		{
			nbbins = (int)(1/x);
			binsname = new String[nbbins];
			double incr = x;
			for(int i = 0; i< binsname.length; i++)
			{
				binsname[i] = ""+round((incr-x),100)+"-"+round(incr,100);
				incr += x;
			}
		}
		else exit("This option (-abs) is not implemented yet!!");
		//{
			//nbbins = max/x;
		//}
		
		
		
		// 2-Read chromosomes files names --> Files + Positions
		
		BufferedReader br;
		String l;
		String[] lT;
		
		ArrayList Files = new ArrayList(); //list of strings
		ArrayList Positions = new ArrayList(); // list of table[int] [0]= centromere [1]= telomere
		
		
		try
		{
			br = new BufferedReader( new FileReader(filelist) );
			l = br.readLine();
			while(l!=null)
			{
				lT= l.split("\t");
				if (lT.length!=3)
				{
					exit("File reading problem : chromosome list must have 3 columns (name \t centromere \t telomere): "+lT.length);
				}
				Files.add(lT[0]);
				int[] tab = new int[2];
				tab[0]= Integer.parseInt(lT[1]);
				tab[1]= Integer.parseInt(lT[2]);
				Positions.add(tab);
				
				l = br.readLine();
			}
			br.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		System.err.println("Reading "+Files.size()+" files ... done");
		
		
		
		// 3-Data Structure

		ArrayList Families = new ArrayList();
		
		//HashMap of String tables
		HashMap MaternalChild = new HashMap(); // [key : name+motherid / values : 22chr+total * nbbins]
		HashMap PaternalChild = new HashMap(); // same but fatherid
		HashMap Parents = new HashMap(); // [key : parent sex nbchild / values : 22chr+total * nbbins] ( or x2 because before and after centromere)
		
		//index
		HashMap IndivToKey =  new HashMap();
		
		//int mapIndex = 0;
		//int mapIndexParents = 0;
		
		// 4-Initialize with famfile
		
		try
		{
			br = new BufferedReader( new FileReader(famfile) );
			l = br.readLine();
			while(l!=null)
			{
				lT = l.split("\t");		//famId \t father \t mother \t child1 child2 child3 ... \t ages
				if(lT.length<4) exit("Pb in famfile");
				Families.add(lT[0]);	//famId
				String[] lTT = lT[3].split(" "); //childrens
		
				int enf = 0;
				for(int i=0; i<lTT.length; i++)
				{
					if( Group==null || (Group!=null && Group.contains(lTT[i])))
					{
					   enf++;
						
					   //childtable are key in hashmap ex: MaternalChild. Here they are initialized 
					   //names + parents
					   String childM = lTT[i]+" "+lT[2];
					   double[][] childtableM = new double[Files.size()+1][nbbins];
					   String childP = lTT[i]+" "+lT[1];
					   double[][] childtableP = new double[Files.size()+1][nbbins];
						
					   //put in child hash map
					   MaternalChild.put(childM,childtableM);
					   PaternalChild.put(childP,childtableP);
						
					   //update index
					   IndivToKey.put(lTT[i], childP+"-"+childM);
					}
				}
				
				if(enf>0)
				{
					//mother
					
					String sM = lT[2]+" female "+ lTT.length; //nb of child
					double[][] tableM = new double[Files.size()+1][nbbins];
					Parents.put(sM,tableM);
					IndivToKey.put(lT[2]+"*",sM);//* to say that it is a parent (to be able to put someone twice if parents and child)
					
					
					//father
					
					String sP = lT[1]+" male "+ lTT.length; //nb of child
					double[][] tableP = new double[Files.size()+1][nbbins];
					Parents.put(sP,tableP);
					IndivToKey.put(lT[1]+"*",sP);//* to say that it is a parent (to be able to put someone twice if parents and child)
				}
				
				/*
				Summary
				
				MaternalChild {child+mother, m}
				PaternalChild {child+father, p}  
				Parents {name+sex+nbchild, par} 
					with m, p and par int[chr+1(tot)][nbbins*2].
				
				IndexToKey {child, child+mother-child+father }
				IndexToKey {mother*, name+sex+nbchild}
				IndexToKey {father*, name+sex+nbchild}
				
				*/
				
				l = br.readLine();
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		if(MaternalChild.size() != PaternalChild.size()) exit("Pb in data structure");
		
		
		
		

		//Chr Data (update using files)
		
		double last = -1;
		
		try
		{
			
			//each file (each chromosome)
			
			for(int chr=0; chr<Files.size(); chr++)
			{
				
				String file = (String)Files.get(chr);
				br = new BufferedReader( new FileReader(file) );
				l = br.readLine();
				
				System.err.print("Reading "+file+ " ...");
				String father = null;
				String mother = null;
				int nbchild = 0;
				
				
				//evaluating lines
				
				while(l!=null)
				{
					//filet
					if(MaternalChild.size() != PaternalChild.size()) exit("Pb in data structure");
					
					lT = l.split(" ");
				
					//Info about the family
					if(lT.length>4 && lT[0].equals("*") && lT[1].equals("Father"))
					{
						father = lT[2];
					}
					else if(lT.length>4 && lT[0].equals("*") && lT[1].equals("Mother"))
					{
						mother = lT[2];
					}
					else if(lT[0].equals("*") && lT[2].equals("children"))
					{
						nbchild = Integer.parseInt(lT[4]);
					}
					
					
					
					//Paternal
					
					else if(lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("paternal"))
					{
						String child = lT[4];
						
						if(IndivToKey.containsKey(child))
						{
							String[] lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in paternal rec");
							
							double tot = Double.parseDouble(lTT[1]); //nb of events
							double[] temp = new double[nbbins];		//events before centromere
							//double[] temp2 = new double[nbbins];	//events after centromere
							l = br.readLine();
							
							if (relative)
							{
								//lines for this child
								for(int i=0;i<tot;i++)
								{
									l = br.readLine();
									lTT = l.split("\t");
									
									int middle = (int)(Double.parseDouble(lTT[2])+((Double.parseDouble(lTT[3])-Double.parseDouble(lTT[2]))/2));
									int[] tabpos = (int[])Positions.get(chr);
									
									//int d = tabpos[0];
									
                                    int index = ComputeRelativeIndex(middle, tabpos[0], tabpos[1], nbbins);
                                    if(index==-1)
                                    {
                                        exit("Position "+middle+" is larger than chromosome last coordinate "+tabpos[1]+"!! Exit.");
                                    }
                                    
									temp[index]++;
                                    
                                    /*if(middle <= tabpos[0] && !signtag.equals("aft")) //before centromere
									{
                                        temp[ComputeRelativeIndex(middle, tabpos[0], d, nbbins)]++;
									}
									else if(middle > tabpos[0] && !signtag.equals("bef"))	//after centromere
									{
										d = tabpos[1]-tabpos[0];
										temp[ComputeRelativeIndex(middle, tabpos[0], d, nbbins)]++;
									}*/
									
									
									//System.err.println(middle +" "+ tabpos[0] +" "+ d +" "+ nbbins);
									//System.err.println(""+ComputeRelativeIndex(middle, tabpos[0], d, nbbins));
								}
							}
							
							
							if(nbchild==2)
							{
								for(int i=0; i<temp.length; i++)
								{
									temp[i] = temp[i]/2; 
									//temp2[i] = temp2[i]/2;
								}
							}
							
							
							//Key for update :
							String fatherkey = null;
							String childkey = null;
							try
							{
								fatherkey = ((String)IndivToKey.get(father+"*"));
								childkey = (((String)IndivToKey.get(child)).split("-"))[0];
							}
							catch(Exception e)
							{
								e.printStackTrace();
								exit("Father ("+father+") not found for child "+child);
							}
							
							//Update :
							double[][] table = (double[][])Parents.get(fatherkey);
							
							for(int n=0; n<table[0].length; n++)
							{
								if(n<temp.length) table[chr][n] += temp[n];	//before centromere
								/*
								else
								{
									table[chr][n] += temp2[n-temp2.length];	//after centromere
								}*/
							}
							
							//exit("");
							
							Parents.remove(fatherkey); //old value
							Parents.put(fatherkey,table); //new value
							
							table = (double[][])PaternalChild.get(childkey);
							table[chr] = temp;
							
							PaternalChild.remove(childkey); //old value
							PaternalChild.put(childkey,table); //new value
						}
						else{
							System.err.println("Child "+child+" not found in fam/group.p");
						}
					}
					
					//Maternal
					
					else if(lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("maternal"))
					{
						String child = lT[4];
						
						if(IndivToKey.containsKey(child))
						{
							String[] lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in maternal rec");
							
							double tot = Double.parseDouble(lTT[1]); //nb of events
							double[] temp = new double[nbbins];
							//double[] temp2 = new double[nbbins];
                            br.readLine();
							
							if (relative)
							{
								last = -1;
								//lines for this child
								for(int i=0;i<tot;i++)
								{
									l = br.readLine();
									lTT = l.split("\t");
									
									if (mb) {
										if(last == -1) last = Double.parseDouble(lTT[3]);
										else
										{
											if(Double.parseDouble(lTT[2])-last < 750000) System.out.println(child+" "+mother+" "+lTT[2]+" "+(chr+1));
											last = Double.parseDouble(lTT[3]);
										}
									}
									
									int middle = (int)(Double.parseDouble(lTT[2])+((Double.parseDouble(lTT[3])-Double.parseDouble(lTT[2]))/2));
									int[] tabpos = (int[])Positions.get(chr);
									
									
									//int d = tabpos[0];
									
									int index = ComputeRelativeIndex(middle, tabpos[0], tabpos[1], nbbins);
                                    if(index==-1)
                                    {
                                        exit("Position "+middle+" is larger than chromosome last coordinate "+tabpos[1]+"!! Exit.");
                                    }
                                    
									temp[index]++;
                                    
                                    /*if(middle <= tabpos[0] && !signtag.equals("aft")) //before centromere
                                     {
                                     temp[ComputeRelativeIndex(middle, tabpos[0], d, nbbins)]++;
                                     }
                                     else if(middle > tabpos[0] && !signtag.equals("bef"))	//after centromere
                                     {
                                     d = tabpos[1]-tabpos[0];
                                     temp[ComputeRelativeIndex(middle, tabpos[0], d, nbbins)]++;
                                     }*/
									
									
									//int d = tabpos[0];
									//if(middle > tabpos[0]) d = tabpos[1]-tabpos[0];
									
									//temp[ComputeRelativeIndex(middle, tabpos[0], d, nbbins)]++;
								}
							}
							
							if(nbchild==2)
							{
								for(int i=0; i<temp.length; i++)
								{
									temp[i] = temp[i]/2; //nb of events
									//temp2[i] = temp2[i]/2; //nb of events
								}
							}
							
							
							//Key for update :
							String motherkey = null;
							String childkey = null;
							try
							{
								motherkey = ((String)IndivToKey.get(mother+"*"));
								childkey = (((String)IndivToKey.get(child)).split("-"))[1];
							}
							catch(Exception e)
							{
								e.printStackTrace();
								exit("Mother ("+mother+") not found for child "+child);
							}
							
							//Update :
							double[][] table = (double[][])Parents.get(motherkey);
							for(int n=0; n<table[0].length; n++)
							{
								if(n<temp.length) table[chr][n] += temp[n]; //before centromere
								/*else
								{
									table[chr][n] += temp2[n-temp2.length]; //after centromere
								}*/
							}
							
							Parents.remove(motherkey); //old value
							Parents.put(motherkey,table); //new value
							
							table = (double[][])MaternalChild.get(childkey);
							//System.err.println(MaternalChild.containsKey(childkey));
							table[chr] = temp;
							
							MaternalChild.remove(childkey); //old value
							MaternalChild.put(childkey,table); //new value
						}
						else{
							System.err.println("Child "+child+" not found in fam/group.m");
						}
					}
					
					l = br.readLine();
				}
				
				System.err.println(" done.");
				br.close();
			}
				
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		if(MaternalChild.size() != PaternalChild.size()) exit("Pb in data structure");
		
		
		
		//update Total
		
		//Maternal + Paternal
		
		Set Keys = MaternalChild.keySet();
		HashMap tempHM1 = new HashMap();
		
		for (Iterator itr = Keys.iterator(); itr.hasNext();)
		//for(int i=0; i<MaternalChild.size(); i++)
		{
			String childkey = (String)itr.next();
			double[][] act = (double[][])MaternalChild.get(childkey);
			
			for(int j = 0; j<act[0].length; j++) //for each bin
			{
				double sum = 0;
				for(int i = 0; i<act.length-1; i++)
				{
					sum += act[i][j];
				}
				act[act.length-1][j] = sum;
			}
			itr.remove();
			tempHM1.put(childkey,act);
			
			/*//replicate
			double[][] act2 = (double[][]) tempHM.get(childkey);
			for(int j = 0; j<act[0].length; j++) //for each bin
			{
				for(int i = 0; i<act.length; i++)
				{
					System.err.print(act[i][j]+" ");
				}
				System.err.println();
			}
			exit("");*/
		}

		MaternalChild = tempHM1;
		//tempHM.clear();
		
		HashMap tempHM2 = new HashMap();
		Keys = PaternalChild.keySet();
		
		for (Iterator itr = Keys.iterator(); itr.hasNext();)
		//for(int i=0; i<MaternalChild.size(); i++)
		{
			String childkey = (String)itr.next();
			double[][] act = (double[][])PaternalChild.get(childkey);
			for(int j = 0; j<act[0].length; j++) //for each bin
			{
				double sum = 0;
				for(int i = 0; i<act.length-1; i++)
				{
					sum += act[i][j];
				}
				act[act.length-1][j] = sum;
			}
			itr.remove();
			
			tempHM2.put(childkey,act);
			
		}
		PaternalChild = tempHM2;
		
		
		//Parents
		
		
		double[] totalP = new double[nbbins]; // total final : to see the proportion of events in each bin.
		double[] totalM = new double[nbbins];
		
		/*
		if (sign) {
			totalP = new double[nbbins*2]; // total final : to see the proportion of events in each bin.
			totalM = new double[nbbins*2];
		}
		*/
		
		HashMap tempHM3 = new HashMap();
		Keys = Parents.keySet();
		
		for (Iterator itr = Keys.iterator(); itr.hasNext();)
		//for(int i=0; i<MaternalChild.size(); i++)
		{
			String namekey = (String)itr.next();
			double[][] act = (double[][])Parents.get(namekey);
			
			for(int j = 0; j<act[0].length; j++) //for each bin
			{
				double sum = 0;
				for(int i = 0; i<act.length-1; i++)
				{
					sum += act[i][j];
					//System.err.print(act[i][j]+" ");
				}
				act[act.length-1][j] = sum;
				//System.err.println(sum);
				
				if(namekey.indexOf("female")!=-1) totalM[j] += sum;
				else totalP[j] += sum;
				/*
				else 
				{
					if (j>=nbbins) 
					{
						if(namekey.indexOf("female")!=-1) totalM[j-nbbins] += sum;
						else totalP[j] += sum;
					}
					else
					{
						if(namekey.indexOf("female")!=-1) totalM[j] += sum;
						else totalP[j] += sum;
					}
				}
				*/

			}
			itr.remove();
			tempHM3.put(namekey,act);
			
		}
		Parents = tempHM3;
		
		
		
		//Printing
		if(stop)
		{
			
			System.out.println("Binnames\tMother\tFather");
			
			double totM = 0;
			double totP = 0;
			
			for(int z=0; z<totalM.length; z++)
			{
				totM += totalM[z];
				totP += totalP[z];
			}
			
			double bins = x;
			for(int z=0; z<totalM.length; z++) 
			{
				//if(Group==null)
				//System.out.println(round(bins,1000)+ "\t"+round(totalM[z]/totM,1000)+"\t"+round(totalP[z]/totP,1000));
				//else 
				
				System.out.println(round(bins,1000)+ "\t"+totalM[z]+"\t"+totalP[z]);
				bins += x;
				
				//System.out.println(binsname[z]+ "\t"+totalM[z]+"\t"+totalP[z]);
			}
			
			
			
			
			exit("Stopped before writing");
		}
		
		//Writing
		
		
		
		BufferedWriter[] BWs = new BufferedWriter[nbbins];
		 
		
		try
		{
			
			
			//.parents FILES
			
			for(int i=0; i<BWs.length; i++)
			{
				String myName = outfile +"_"+binsname[i]+signtag+".parents";
				BWs[i] = new BufferedWriter( new FileWriter( myName ) );
				BWs[i].write("parent\tsex\tnbchildren\t");
				//legend (with chr file names)
				for(int z=0; z<Files.size(); z++)
				{
					String[] spl = ((String)Files.get(z)).split("/");
					BWs[i].write(spl[spl.length-1]+"\t");
				}
				BWs[i].write("total\n");
				//BWs[0].write("allo\n");
				//BWs[0].close();
			}
			
			//for(int i=0; i<BWs.length; i++)
			//{
				//BWs[i].write("zinzan\n");
				//System.err.println(binsname[0]);
			//}
			
			//for(int i=0; i<BWs.length; i++)
			//{
				//BWs[i].close();
				//System.err.println(binsname[0]);
			//}
			//variables to compute moyenne parentale
			HashMap Par = new HashMap();
			String name = "";
			double nb = 0;
			
			
			Keys = Parents.keySet();
			//System.out.println(Parents.size());
			for (Iterator itr = Keys.iterator(); itr.hasNext();)
			{
				//BWs[0].write("test close\n");
				//System.out.println("allo");
				String namekey = (String)itr.next();
			
				double[][] act = (double[][])Parents.get(namekey);
				
				if (act[0].length != BWs.length) exit("Pb in data structure : BWs vs Data");
				
				for(int j = 0; j<act[0].length; j++)	// a travers les diff BWs/bins
				{
					
					name = (namekey.split(" "))[0];
					nb = Double.parseDouble((namekey.split(" "))[2]);
					//if(name.equals("44_8")) exit("ouii");
					
					for(int i = 0; i<(namekey.split(" ")).length; i++)		
					{
						BWs[j].write((namekey.split(" "))[i]+"\t");
					}
					
					for(int i = 0; i<act.length; i++)		// a travers les chr
					{
						BWs[j].write(act[i][j]+"\t");
						
						//computation of moyenne parentale
						if(i==act.length-1) 
						{
							double temp = (double)act[i][j] / nb;
							if(Par.containsKey(name))
							{
								double[] a = (double[])Par.get(name);
								a[j] = temp;
								Par.put(name,a);
							}
							else
							{
								double[] a = new double[nbbins];
								a[j] = temp;
								Par.put(name,a);
							}
						}
					}
					BWs[j].write("\n");
				
				}
			}
			
			for(int i=0; i<BWs.length; i++)
			{
				BWs[i].close();
			}
			
			
			//.females FILES
		
		
			for(int i=0; i<BWs.length; i++)
			{
				String myName = outfile +"_"+binsname[i]+signtag+".females";
				BWs[i] = new BufferedWriter( new FileWriter( myName ) );
				BWs[i].write("child\tparent\t");
				//legend (with chr file names)
				for(int z=0; z<Files.size(); z++)
				{
					String[] spl = ((String)Files.get(z)).split("/");
					BWs[i].write(spl[spl.length-1]+"\t");
				}
				BWs[i].write("total\tpar\n");
			}
			
			
			Keys = MaternalChild.keySet();
			
			for (Iterator itr = Keys.iterator(); itr.hasNext();) // a travers les individus
			{
				String namekey = (String)itr.next();
				double[][] act = (double[][])MaternalChild.get(namekey);
				name = (namekey.split(" "))[1];
				double[] a = (double[])Par.get(name); //parental means for each bins
				
				if (act[0].length != BWs.length) exit("Pb in data structure : BWs vs Data");
				
				for(int j = 0; j<act[0].length; j++)	// a travers les diff BWs/bins
				{
					
					for(int i = 0; i<(namekey.split(" ")).length; i++)		
					{
						BWs[j].write((namekey.split(" "))[i]+"\t");
					}
					
					for(int i = 0; i<act.length; i++)		// a travers les chr
					{
						BWs[j].write(act[i][j]+"\t");
					}
					
					BWs[j].write(round(a[j],1000)+"\n");
				}
			}
			
			for(int i=0; i<BWs.length; i++)
			{
				BWs[i].close();
			}
			
			

			
			
			//.males FILES
		
		
			for(int i=0; i<BWs.length; i++)
			{
				String myName = outfile +"_"+binsname[i]+signtag+".males";
				BWs[i] = new BufferedWriter( new FileWriter( myName ) );
				BWs[i].write("child\tparent\t");
				//legend (with chr file names)
				for(int z=0; z<Files.size(); z++)
				{
					String[] spl = ((String)Files.get(z)).split("/");
					BWs[i].write(spl[spl.length-1]+"\t");
				}
				BWs[i].write("total\tpar\n");
			}
			
			
			Keys = PaternalChild.keySet();
			
			for (Iterator itr = Keys.iterator(); itr.hasNext();) // a travers les individus
			{
				String namekey = (String)itr.next();
				double[][] act = (double[][])PaternalChild.get(namekey);
				name = (namekey.split(" "))[1];
				
				double[] a = (double[])Par.get(name); //parental means for each bins
				if(a == null) exit("pb in data structure : "+namekey);
				
				if (act[0].length != BWs.length) exit("Pb in data structure : BWs vs Data");
				
				for(int j = 0; j<act[0].length; j++)	// a travers les diff BWs/bins
				{
					
					
					for(int i = 0; i<(namekey.split(" ")).length; i++)		
					{
						BWs[j].write((namekey.split(" "))[i]+"\t");
					}
					
					for(int i = 0; i<act.length; i++)		// a travers les chr
					{
						BWs[j].write(act[i][j]+"\t");
					}
					
					BWs[j].write(round(a[j],1000)+"\n");
				}
			}
			
			for(int i=0; i<BWs.length; i++)
			{
				BWs[i].close();
			}
			
			

		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		
		
	}

	//Utilitaires
	
	public static void exit(String s)
	{
	
		System.err.println(s);
		System.exit(0);
	
	}
	
	public static double round(double a, int level)
	{
		int a2 = (int)(a * level + 0.5);
		
		return (double)a2/level;
		
		
	}
	
	public static double ComputeMean(double[] table)
	{
		double sum = 0;	
		for(int i=0; i<table.length; i++)
		{
			sum += table[i];
		}
		//System.err.println(sum+" "+table.length);
		return sum/table.length;
	}
	
	//return index if every things ok
	//return -1 if val > d (makes no sense)
	//return -2 if there is a pb in the returning of the index part
    
    public static int ComputeRelativeIndex(int val, int c, int d, int nbbins)
	{
		double x = 0; //relative value
		
		int d1 = c;
		int d2 = d-c;
		
		if (val > d) return -1;
		if(val < c) // the event is btw 0 and centromere
		{
			x = (double)(d1-val)/d1;
		}
		else
		{
			x = (double)(d2-(d-val))/d2;
		}
		//System.err.println(x);
		
		
		//small bins are closer to centromere
		double binsize = 1/(double)nbbins;
		
		double incr = binsize;
		for(int i=0; i<nbbins; i++)
		{
			if(x<=incr) return i;
			else incr += binsize;
			//System.err.println(incr);
		}
		
		return -2;
	}
    
    //earlier version this one is wrong
	public static int ComputeRelativeIndexWrong(int val, int c, int d, int nbbins)
	{
		double x = 0; //relative value
		
		if(val < c) // the event is btw 0 and centromere
		{
			if (val > d) x = -1;
			else x = (double)(d-val)/d;
		}
		else
		{
			x = (double)(val-c)/d;
		}
		//System.err.println(x);
		
		double binsize = 1/(double)nbbins;
		
		double incr = binsize;
		for(int i=0; i<nbbins; i++)
		{
			if(x<incr) return i;
			else incr += binsize;
			//System.err.println(incr);
		}
		
		return -2;
	}
	

}