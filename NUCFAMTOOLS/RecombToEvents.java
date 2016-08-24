//
//  RecombToEvents.java
//  
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Summarize recombination events form Recomb to a format for Overlap program
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class RecombToEvents {
	
	
	public static void main(String[] args)
	{
		
		
		
		if (args.length < 2) exit(" usage : java -cp NucFamTools.jar RecombToEvents -l/-f file.recomb [-o outfile -fam excludeFamilies.txt -i individuals.txt -chrfile -printk]");
		
		//list in order of chr (in the file.events, the rank of the file will be in the chr column) 
		//EXCEPT if -chrfile option is specified :
		//in that case, the chr index will be taken from file name : the chr index should be preceded by chr (and should be followed by "_")
		
		
		
		
		
		ArrayList Files = new ArrayList();
		if(args[0].equals("-f"))
		{
			Files.add(args[1]);
		}
		else if(args[0].equals("-l"))
		{
			try
			{
				BufferedReader br = new BufferedReader( new FileReader(args[1]) );
				String l = br.readLine();
				while(l!=null)
				{
					Files.add(l);
					l = br.readLine();
				}
				br.close();
			}
			catch(Exception e)
			{
				e.printStackTrace();
				exit("File reading problem");
			}
		}
		else
		{
			exit(" usage : java -cp NucFamTools.jar RecombToOverlap -l/-f file.recomb [-o outfile -fam excludeFamilies.txt -i individuals.txt -chrfile -printk]");
		}
		
		System.err.println("Charging "+Files.size()+" files ... done");
		
		
		TreeSet Families = new TreeSet();		//families to exclude
		TreeSet Indiv = new TreeSet();			//individuals to include in remaining families (to evaluate maternal age)
		String outfile = args[1]+".events";
		
		boolean newNFrecomb = true;			//changes in april 2011 for cancer project, format of output changed
		boolean chrfile = false;
		boolean printk = false;
		
		BufferedReader br;
		String l;
		String[] lT;
		
		BufferedWriter bw;
		
		if(args.length > 2)
		{
			int index = 2;
			while(index < args.length)
			{
				if(args[index].equals("-o"))
				{
					outfile = args[index+1];
					index += 2;
				}
				else if(args[index].equals("-chrfile")) //option added 02/09/2011
				{
					chrfile = true;
					index += 1;
				}
				else if(args[index].equals("-printk")) //option added 09/09/2011
				{
					printk = true;
					index += 1;
				}
				else if(args[index].equals("-fam"))
				{
					try
					{
						br = new BufferedReader( new FileReader(args[index+1]) );
						l = br.readLine();
						while(l!=null)
						{
							Families.add(l);
							l = br.readLine();
						}
						index += 2;
					}
					catch(Exception e)
					{
						e.printStackTrace();
						exit("File reading problem");
					}
				}
				else if(args[index].equals("-i"))
				{
					try
					{
						br = new BufferedReader( new FileReader(args[index+1]) );
						l = br.readLine();
						while(l!=null)
						{
							Indiv.add(l);
							l = br.readLine();
						}
						index += 2;
					}
					catch(Exception e)
					{
						e.printStackTrace();
						exit("File reading problem");
					}
				}
				else
				{
					exit("Invalid option : "+ args[index]);
					
				}
			}
			
			
		}
		
		
		
		try
		{
			
			bw = new BufferedWriter( new FileWriter( outfile ) );
			bw.write("chr\tleft\tright\tsex\tsmall_less_than_4\twidth\tchild\tfam");
			if(printk)bw.write("\tk_bef\tk_aft");
			bw.newLine();
			
			
			for(int chr=1; chr<=Files.size(); chr++)
			{
				String file = (String)Files.get(chr-1);
				String id = ""+chr;
				
				//System.err.println("allo");
				if(chrfile){ id = chrid(file);}
				
				
				br = new BufferedReader( new FileReader(file) );
				l = br.readLine();
				
				System.err.print("Reading "+file+ " ...");
				int nbchild = -1;
				boolean lookForF = true;
				boolean lookForM = true;
				
				boolean thisFam = true;
				String fam="";
				while(l!=null)
				{
					
					lT = l.split(" ");
					
					if(lT[0].equals("LOADING"))
					{
						fam = lT[2];
						if (Families.contains(fam))
						{
							thisFam = false;
						}
						else thisFam = true;
					}
					
					else if(thisFam && lT.length>4 && lT[0].equals("*") && lT[2].equals("children"))
					{
						nbchild = Integer.parseInt(lT[4]);
						lookForF = true;
						lookForM = true;
					}
					
					
					
					//meanEventF
					else if(thisFam && lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("paternal"))
					{
						
						if ((lookForF && Indiv.contains(lT[4])) || (lookForF && Indiv.size()==0) )
						{
							//System.err.println("z");
							String[] lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in meanEventF");
							double temp = Double.parseDouble(lTT[1]); //nb of events
							
							if(newNFrecomb) l = br.readLine();
							//lines
							for(int i=0;i<temp;i++)
							{
								l = br.readLine();
								lTT = l.split("\t");
								
								int diff = (Integer.parseInt(lTT[3])-Integer.parseInt(lTT[2]));
								String small = "small_family";
								if(nbchild > 3) small = "not_small_family";
								bw.write(id+"\t"+lTT[2]+"\t"+lTT[3]+"\tmale\t"+small+"\t"+diff+"\t"+lT[4]+"\t"+fam);//+fam added 8/10/2011
								if(printk)bw.write("\t"+lTT[1]+"\t"+lTT[4]);
								bw.newLine();
								
							}
							
							
							if(nbchild == 2)
							{
								lookForF = false;
							}
						}
					}
					
					//meanEventM
					else if(thisFam && lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("maternal"))
					{
						if ((lookForM && Indiv.contains(lT[4])) || (lookForM && Indiv.size()==0) )
						{
							String[] lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in meanEventM");
							double temp = Double.parseDouble(lTT[1]);
							
							if(newNFrecomb) l = br.readLine();
							//lines
							for(int i=0;i<temp;i++)
							{
								l = br.readLine();
								lTT = l.split("\t");
								
								int diff = (Integer.parseInt(lTT[3])-Integer.parseInt(lTT[2]));
								String small = "small_family";
								if(nbchild > 3) small = "not_small_family";
								bw.write(id+"\t"+lTT[2]+"\t"+lTT[3]+"\tfemale\t"+small+"\t"+diff+"\t"+lT[4]+"\t"+fam);//+fam added 8/10/2011
								if(printk)bw.write("\t"+lTT[1]+"\t"+lTT[4]);
								bw.newLine();
								
							}
							
							
							if(nbchild == 2)
							{
								lookForM = false;
							}
						}
					}
					
					
					l = br.readLine();
				}
				
				System.err.println(" done.");
				br.close();
			}
			
			bw.close();	
			System.err.println("Writing "+outfile+" ... done.");		
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
	
	public static String chrid(String s)
	{
		int indexc = s.indexOf("chr");
		String s2 = s.substring(indexc+3,indexc+5);
		s2 = s2.replaceAll("_","");
		return s2;
	}
	
}