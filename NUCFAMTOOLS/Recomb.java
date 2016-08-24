//
//  Recomb.java
//  (initially : NF_Recomb.java in LVTOOLS)
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Main Program (using NF classes)
//
//	throughout the program F is used for Father (male)
//						   M is used for Mother (female)

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class Recomb {


	public static TreeSet ValidErrors;
	public static boolean cbool;
	public static boolean sbool;
	
	public static void main(String[] args)
	{
	
		int minArg = 2;
		if (args.length < minArg) exit(" usage : java -cp NucFamTools.jar Recomb file.ped file.map [-o outfile -x -k N -err value(Kb) -valerr file.err -sharing -nocrossover]");
		
		String pedfile = args[0];
		String mapfile = args[1];
		String outfile = null;
		int Kprint = 2;
		int interval = 10000000;//pb
		ArrayList Families = new ArrayList();
		boolean chrx = false;
		ValidErrors = null;
		
		cbool = true;
		sbool = false;
		
		if(args.length > minArg)
		{
			int index = minArg;
			while(index < args.length)
			{
				if(args[index].equals("-err"))
				{
					interval = Integer.parseInt(args[index+1]) * 1000;//pb
					index += 2;
				}
				else if(args[index].equals("-k"))
				{
					Kprint = Integer.parseInt(args[index+1]);
					index += 2;
				}
				else if(args[index].equals("-x"))
				{
					chrx = true;
					index += 1;
				}
				else if(args[index].equals("-nocrossover"))
				{
					cbool = false;
					index += 1;
				}
				else if(args[index].equals("-sharing"))
				{
					sbool = true;
					index += 1;
				}
				else if(args[index].equals("-o"))
				{
					outfile = args[index+1]+"_";
					index += 2;
				}
				else if(args[index].equals("-valerr"))
				{
					try
					{
						BufferedReader br = new BufferedReader( new FileReader(args[index+1]) );
						ValidErrors = new TreeSet();
						String l = br.readLine();
						while(l!=null)
						{
							ValidErrors.add(l);
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
		if (outfile == null)
		{
			if(interval!=10000000)
			{
				outfile = pedfile+"_"+(interval/1000)+"kb_";
			}
			else
			{
				outfile = pedfile+"_";
			}
		}
		
		if(!cbool && !sbool)
		{
			exit("-nocrossover option is specified without the -sharing option. There is nothing left to do here! Exit.");
		
		}
		
		if(sbool)
		{
			File f = new File(outfile+"sharing.txt");
			if(f.exists()) 
			{ 
				try 
				{
					BufferedWriter bw = new BufferedWriter( new FileWriter( outfile+"sharing.txt" ) );
					//just creating it
					bw.close();
				}
				catch(Exception e)
				{
					e.printStackTrace();
					exit("File writing problem (sbool).");
				}
			}
		}
		
		
	
		try 
		{
			BufferedWriter bw = new BufferedWriter( new FileWriter( outfile+"infomarkF.txt" ) );
			bw.write("MarkerNbF\tPosition\tFam\tChildrensAlleles\n");
			bw.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File writing problem (infomarkersF).");
		}

		
		try 
		{
			BufferedWriter bw = new BufferedWriter( new FileWriter( outfile+"infomarkM.txt" ) );
			bw.write("MarkerNbM\tPosition\tFam\tChildrensAlleles\n");
			bw.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File writing problem (infomarkersM).");
		}
		
		
		//System.out.println("java -cp NucFamTools.jar Recomb "+args[0]+" "+args[1]+" -err "+interval+" -k "+ Kprint);
		
		
		
		BufferedReader br;
		String l;
		String[] lT;
		
		
		
		
		try
		{
			
			br = new BufferedReader( new FileReader(mapfile) );
			ArrayList Positions = new ArrayList();
			
			l = br.readLine();
			while(l!=null)
			{
				Positions.add(l);
				l = br.readLine();
			}
			
			
			br = new BufferedReader( new FileReader(pedfile) );
			NF_NucFam Family;
			ArrayList maFamille;
			
			l = br.readLine();	//first line of first family
			boolean end = false;
			
			while(!end)
			{
				String id;
				int nbchildren;
				maFamille = new ArrayList();
				
				lT = l.split("\t");
				if(lT.length!=2 && !lT[0].equals("Family ID :")) exit("Problem in Family ID line in file "+pedfile+".");
				id = lT[1];
				
				
				l = br.readLine();
				lT = l.split("\t");
				if(lT.length!=2 && !lT[0].equals("Number of children :")) exit("Problem in Number of children line in file "+pedfile+".");
				nbchildren = Integer.parseInt(lT[1]);
				
				l = br.readLine(); //empty line
				
				//Data
				for(int i=0; i<nbchildren+2; i++) //family = 2 parents + childrens
				{
					l = br.readLine();
					lT = l.split("\t");
					if(lT.length<7)	exit("Problem in Data Format for file "+pedfile+" (less than 7 obligatory columns).");
					if(lT[6].length() !=3) 	exit("Problem in Data Format for file "+pedfile+" (wrong format for first locus : "+lT[6]+").");
					
					maFamille.add(l);
				}
				
				//Family
				
				Family = new NF_NucFam(id, maFamille, Positions, interval, Kprint, chrx, outfile);
				//What does it do : 
				// - verify family format (father/mother/children) and loci format
				// - SetHapInformative() for parents
				// - children in TreeSet
				// - Recode children (class NF_Child)
				// - Find errors
				// - Crossovers
				Families.add(Family.Printing);
				
				
				l = br.readLine(); //empty line
				if(l==null) end = true;
				else
				{
					l = br.readLine(); //first line of next family
					if(l==null) end = true;
				}
			}
			
			br.close();
			
			if(cbool)
			{
			
			for(int i=1; i<Kprint+2; i++)
			{
				BufferedWriter bw = new BufferedWriter( new FileWriter( outfile+"k"+(i-1)+".recomb" ) );
				for(int f = 0;f<Families.size(); f++)
				{
					bw.write("LOADING FAMILY "+(f+1)+" ... \n");
					String[] prt = (String[])Families.get(f);
					bw.write(prt[0]);
					bw.write(prt[i]);
				} 
				bw.close();
				
			}
				
			}
			
			if(sbool)
			{
				System.err.println("Writing sharing between siblings in "+outfile+ "sharing.txt ... done");
			}
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File reading problem");
		}
		
		System.err.println("Bye bye!\n");
	
	}

	//Utilitaires
	
	public static void exit(String s)
	{
	
		System.err.println(s);
		System.exit(0);
	
	}
	

}