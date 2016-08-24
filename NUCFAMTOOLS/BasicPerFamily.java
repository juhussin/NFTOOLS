//
//  BasicPerFamily.java
//  NUCFAMTOOLS
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	First step analysis : recap : informative markers, rec events, errors --> per family.
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class BasicPerFamily {


	public static void main(String[] args)
	{
		if (args.length < 3) exit(" usage : java -cp NucFamTools.jar BasicPerFamily nbfamily -l/-f file.recomb [-o outfile.res]");
		
		/*System.err.println(args.length);
		for(int i =0; i<args.length;i++)
		{
			System.err.println(i+" "+args[i]);
		}*/
		
		int nbFam = Integer.parseInt(args[0]);
		boolean list = true;
		if(args[1].equals("-f"))
		{
			list = false;
		}
		else if(args[1].equals("-l"))
		{
			list = true;
		}
		else
		{
			exit(args[1]+ "is not a valid file option (must be -l for filelist or -f for a unique file");
		}
		
		String file = args[2];
		String outfile = args[2]+".res";
		
		if(args.length > 3 && args[3].equals("-o"))
		{
			outfile = args[4]; 
		}

		
		ArrayList Files =  new ArrayList();
		BufferedReader br;
		String l;
		String[] lT;

		BufferedWriter bw;
		
		
		if(list)
		{
			try
			{
				br = new BufferedReader( new FileReader(file) );
				l = br.readLine();
				while(l!=null)
				{
					Files.add(l);
					l = br.readLine();
				}
			}
			catch(Exception e)
			{
				e.printStackTrace();
				exit("File reading problem");
			}
		}
		else
		{
			Files.add(file);
		}
		
		
		
		String[][] ResName = new String[nbFam][4];
		boolean donename = false;
		int[][] ResData = new int[nbFam][6];
		

		for(int ia = 0; ia<Files.size(); ia++)
		{
		
		
		try
		{
			br = new BufferedReader( new FileReader((String)Files.get(ia)) );
			//l = br.readLine();
			
			System.err.print("Reading "+(String)Files.get(ia)+ " ...");
			int nbchild = -1;
			//boolean lookFor = true;
			
			l = br.readLine();
			lT = l.split(" ");
			
			for(int fam=0; fam<nbFam; fam++)
			{
				if (!lT[0].equals("LOADING")) exit("pb 0");
				if (ResName[fam][0]==null) ResName[fam][0] = lT[2];
				
				//father
				while(!(lT.length>4 && lT[0].equals("*") && lT[1].equals("Father"))){
					l= br.readLine();
					lT = l.split(" ");
				}
				ResData[fam][1] += Integer.parseInt(lT[4]);
				if (ResName[fam][1]==null) ResName[fam][1] = lT[2];
				
				//mother
				l= br.readLine();
				lT = l.split(" ");
				if(!(lT.length>4 && lT[0].equals("*") && lT[1].equals("Mother")))
				{
					exit("pb 1");
				}
				else
				{
					ResData[fam][2] += Integer.parseInt(lT[4]);
					if (ResName[fam][2]==null) ResName[fam][2] = lT[2];
				}
				
				//children
				l= br.readLine();
				lT = l.split(" ");
				if(!(lT.length>4 && lT[0].equals("*") && lT[2].equals("children")))
				{
					exit("pb 2");
				}
				else
				{
					if (ResData[fam][0]==0) ResData[fam][0] = Integer.parseInt(lT[4]);
				}
				
				l = br.readLine();
				lT = l.split("\t");
				
				//errors
				while(!(lT.length>1 && lT[1].equals("Total:"))){ //changed 8/10/2011 lT[0].equals("TOTAL:") --> lT[1].equals("Total:")
					l= br.readLine();
					lT = l.split("\t");
					//System.err.println(lT[0]);
				}
				String[] lTT = lT[2].split(" "); //changed 8/10/2011 lT[1].split(" ")
				int nberror = Integer.parseInt(lTT[0]);
				ResData[fam][5] += nberror;
				
				l = br.readLine();
				lT = l.split(" ");
				
				//rec events
				boolean lookforF = true;
				boolean lookforM = true;
				
				while(true)
				{
					l= br.readLine();
					if(l==null) break;
					
					lT = l.split(" ");
					if(lT[0].equals("LOADING")) break;
					
					if(lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("paternal"))
					{
						if(!donename)
						{
							if(ResName[fam][3]==null) ResName[fam][3] = lT[4];
							else ResName[fam][3] = ResName[fam][3]+"/"+lT[4];
						}
						if (lookforF)
						{
							lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in meanEventF");
							int temp = Integer.parseInt(lTT[1]);
							
							ResData[fam][3] += temp;
							
							if(ResData[fam][0] == 2)
							{
								lookforF = false;
							}
						}
					}
					
					//meanEventM
					if(lT.length> 10 && lT[0].equals("Crossover") && lT[9].equals("maternal"))
					{
						if (lookforM)
						{
							lTT = lT[lT.length-1].split("\t");
							if(lTT.length!=2) exit("pb in meanEventM");
							int temp = Integer.parseInt(lTT[1]);
							
							ResData[fam][4] += temp;
							
							if(ResData[fam][0] == 2)
							{
								lookforM = false;
							}
						}
					}
				}
			}
			
			
			
			System.err.println(" done.");
			br.close();
			donename =true;			
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		}
		
		try
		{
			bw = new BufferedWriter( new FileWriter(outfile) );
			bw.write("#Summary per family\n#\nFamily\tFather\tMother\tChildrensNames\tnbchild\tIMfather\tIMmother\tRECpaternal\tRECmaternal\tErrors\tRECpatPerChild\tRECmatPerChild\n");
			
			double meanF = 0;
			double meanM = 0;
			double infoF = 0;
			double infoM = 0;
			double meiosis = 0;
			
			for(int i=0; i<ResData.length; i++)
			{
				for(int j=0; j<ResName[0].length; j++)
				{
					bw.write(ResName[i][j]+"\t");
				}
				for(int j=0; j<ResData[0].length; j++)
				{
					bw.write(ResData[i][j]+"\t");
				}
				bw.write(round(ResData[i][3]/(double)ResData[i][0],100)+"\t"+round(ResData[i][4]/(double)ResData[i][0],100)+"\n");
				
				meiosis += ResData[i][0];
				infoF += ResData[i][1];
				infoM += ResData[i][2];
				meanF += ResData[i][3];
				meanM += ResData[i][4];
			}
			
			System.err.println( "Number of meiosis : "+ meiosis);
			System.err.println( "Mean number of Paternal informative markers : "+ round(infoF/ResData.length,1000));
			System.err.println( "Mean number of Maternal informative markers : "+ round(infoM/ResData.length,1000));
			System.err.println( "Mean number of Paternal recombination per meiosis : "+ round(meanF/meiosis,1000));
			System.err.println( "Mean number of Maternal recombination per meiosis : "+ round(meanM/meiosis,1000));
			
			System.err.println();
			System.err.println("Writing "+outfile+" ... done");
			bw.close();
			
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File writing problem");
		}
		
		

		
	
	}

	//Utilitaires
	
	public static void exit(String s)
	{
	
		System.err.println(s);
		System.exit(0);
	
	}
	
	//level = 1, 10, 100, 1000, 10000
	public static double round(double a, int level)
	{
		int a2 = (int)(a * level + 0.5);
		
		return (double)a2/level;
		
	}

}