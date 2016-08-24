//
//  (Initially : NF_LinkageToRecomb.java in LVTOOLS)
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Conversion Program : Takes Linkage files and convert it in a format understood by NF_Recomb
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class PlinkToRecomb {
	
	
	public static void main(String[] args)
	{
		/*
		
		if (args.length < 1) exit(" usage : java -cp NucFamTools.jar PlinkToRecomb file.ped [-o outfile -orderedbyfam ]");
		
		String pedfile = args[0];
		String outfile = args[0].replaceAll(".ped","_recomb.ped");
		
		if(args.length > 1)
		{
			outfile = args[2]; 
		}
		else if (outfile.equals(pedfile)) //no .ped in the pedfile name
		{
			outfile = args[0]+"_recomb.ped";
		}
		*/
		
		int minArg = 1;
		if (args.length < minArg) exit(" usage : java -cp NucFamTools.jar PlinkToRecomb file.ped [-o outfile -nbid]");
		
		String pedfile = args[0];
		
		String outfile = args[0].replaceAll(".ped","_recomb.ped");
		
		boolean orderedbyfam = false;
		boolean nbid = false;
		
		
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
				else if(args[index].equals("-ordByFam"))
				{
					orderedbyfam = true;
					index += 2;
				}
				else if(args[index].equals("-nbid"))
				{
					nbid = true;
					index += 2;
				}
				else
				{
					exit("Invalid option : "+ args[index]);
					
				}
			}
			
			
		}
		
		if (outfile.equals(pedfile)) //no .ped in the pedfile name and no -o option
		{
			outfile = args[0]+"_recomb.ped";
		}
		
		
		
		BufferedReader br;
		String l;
		String[] lT;
		
		BufferedWriter bw;
		
		TreeMap ParentsToChildren = new TreeMap();
		TreeMap Data = new TreeMap();
		
		try
		{
			br = new BufferedReader( new FileReader(pedfile) );
			l = br.readLine();
			int cc = 0;
			
			System.err.print("Reading "+pedfile+ " ");
			
			while(l!=null)
			{
				if(cc%50 == 0) System.err.print(".");
				
				
				lT = l.split("\t");
				if(lT.length<7)	exit("Problem in Data Format for file "+pedfile+" (there are less than 7 obligatory columns).");
				if(lT[6].length() !=3) 	exit("Problem in Data Format for file "+pedfile+" (wrong format for first locus : "+lT[6]+").");
				
				//Data
				
				//changes 13/11/2013
				//not reading the entire file anymore, juste keeping the line : cc
				/*
				String line = lT[1]+"-"+lT[4]+"-";
				for(int i=6; i<lT.length-1; i++)
				{
					line += lT[i]+"\t"; 
				}
				line += lT[lT.length-1];
				 
				Data.put(lT[1], line);
				*/
				
				
				
				Data.put(lT[1], new Integer(cc));
				
				if(!lT[2].equals("0") && !lT[3].equals("0"))
				{
					String nucfam = lT[2]+"-"+lT[3];
					if(ParentsToChildren.containsKey(nucfam))
					{
						String act = (String)ParentsToChildren.get(nucfam);
						ParentsToChildren.put(nucfam, act+"-"+lT[1]);
					}
					else
					{
						ParentsToChildren.put(nucfam, lT[1]);
					}
				}
				
				l = br.readLine();
				cc++;
			}
			
			System.err.println(" done.");
			br.close();
			
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		
		System.err.print("Compiling nuclear families with 2 and more children ... ");
		
		Set Parents = ParentsToChildren.keySet();
		TreeMap FatherToMother = new TreeMap();
		for (Iterator itr = Parents.iterator(); itr.hasNext();)
		{
			String act = (String)itr.next();
			String[] acT = act.split("-");
			if(Data.containsKey(acT[0]) && Data.containsKey(acT[1]))
			{
				String[] children = ((String)ParentsToChildren.get(act)).split("-");
				if(children.length > 1)
				{
					FatherToMother.put(acT[0], acT[1]);
				}
			}
		}
		
		System.err.println(" done ("+FatherToMother.size()+" nuclear families found).");
		
		
		try
		{
			
			System.err.print("Writing "+outfile+ " ");
			
			
			bw = new BufferedWriter( new FileWriter( outfile ) );
			int countfam = 1;
			
			Set Fathers = FatherToMother.keySet();
			
			int cc =0;
			
			for (Iterator itr = Fathers.iterator(); itr.hasNext();)
			{
				if(cc%10 == 0) System.err.print(".");
				cc++;
				
				
				String father = (String)itr.next();
				String mother = (String)FatherToMother.get(father);
				String parents = father+"-"+mother;
				
				String childrenline = (String)ParentsToChildren.get(father+"-"+mother);
				
				if(childrenline == null)
				{
					exit("Map problem (should never happen)");
				}
				
				String[] children = childrenline.split("-");
				
				if(nbid)bw.write("Family ID :\t"+countfam+"\n");
				else bw.write("Family ID :\t"+parents+"\n");
				bw.write("Number of children :\t"+children.length+"\n\n");
				//changes 13/11/2013
				//not reading the entire file anymore, juste keeping the line : cc
				//needs to fetch the line in pedfile again
				/*
				bw.write((String)Data.get(father)+"\n");
				bw.write((String)Data.get(mother)+"\n");
				for(int i=0; i<children.length; i++)
				{
					bw.write((String)Data.get(children[i])+"\n");
				}
				*/
				
				//father
				int index = ((Integer)Data.get(father)).intValue();
				BufferedReader brtemp = new BufferedReader(new FileReader(pedfile));
				String ltemp = brtemp.readLine();
				
				for(int i=0;i<index+1;i++)
				{
					if(i==index)
					{
						lT = ltemp.split("\t");
						
						bw.write(lT[1]+"-"+lT[4]+"-");
						for(int it=6; it<lT.length-1; it++)
						{
							bw.write(lT[it]+"\t"); 
						}
						bw.write(lT[lT.length-1]+"\n");
						
					}
					ltemp = brtemp.readLine();
				}
				
				//mother
				int previndex = index;
				index = ((Integer)Data.get(mother)).intValue();
				
				if(previndex<index)
				{
					for(int i=previndex+1;i<index+1;i++)
					{
						if(i==index)
						{
							lT = ltemp.split("\t");
							
							bw.write(lT[1]+"-"+lT[4]+"-");
							for(int it=6; it<lT.length-1; it++)
							{
								bw.write(lT[it]+"\t"); 
							}
							bw.write(lT[lT.length-1]+"\n");
							
						}
						ltemp = brtemp.readLine();
					}
				}
				else
				{
					//start over
					brtemp.close();
					brtemp = new BufferedReader(new FileReader(pedfile));
					ltemp = brtemp.readLine();
					
					for(int i=0;i<index+1;i++)
					{
						if(i==index)
						{
							lT = ltemp.split("\t");
							
							bw.write(lT[1]+"-"+lT[4]+"-");
							for(int it=6; it<lT.length-1; it++)
							{
								bw.write(lT[it]+"\t"); 
							}
							bw.write(lT[lT.length-1]+"\n");
							
						}
						ltemp = brtemp.readLine();
					}
				
				}
				
				//children
				for(int ch=0; ch<children.length; ch++)
				{
					previndex = index;
					index = ((Integer)Data.get(children[ch])).intValue();
					
					//same as mother 
					
					if(previndex<index)
					{
						for(int i=previndex+1;i<index+1;i++)
						{
							if(i==index)
							{
								lT = ltemp.split("\t");
								
								bw.write(lT[1]+"-"+lT[4]+"-");
								for(int it=6; it<lT.length-1; it++)
								{
									bw.write(lT[it]+"\t"); 
								}
								bw.write(lT[lT.length-1]+"\n");
								
							}
							ltemp = brtemp.readLine();
						}
					}
					else
					{
						//start over
						brtemp.close();
						brtemp = new BufferedReader(new FileReader(pedfile));
						ltemp = brtemp.readLine();
						
						for(int i=0;i<index+1;i++)
						{
							if(i==index)
							{
								lT = ltemp.split("\t");
								
								bw.write(lT[1]+"-"+lT[4]+"-");
								for(int it=6; it<lT.length-1; it++)
								{
									bw.write(lT[it]+"\t"); 
								}
								bw.write(lT[lT.length-1]+"\n");
								
							}
							ltemp = brtemp.readLine();
						}
						
					}
					
					
				}
				
				brtemp.close();
				
				bw.write("\n");
				countfam++;
			}
			
			Data.clear();
			Fathers.clear();
			FatherToMother.clear();
			ParentsToChildren.clear();
			
			System.err.println(" done.");
			
			bw.close();
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
	
	
}