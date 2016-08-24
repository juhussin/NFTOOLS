//
//  NF_NucFam.java
//  LVTOOLS
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Object NF_NucFam
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class NF_NucFam {
	
	//SIMPLE NUCLEAR FAMILY
	
	
	
	//RECOMB
	
	String ID;
	NF_Parent Father;
	NF_Parent Mother;
	
	TreeSet Children;
	
	int[] AllPositions;
	TreeSet AllErrors;
	String[] Printing;
	int locisize = 0;
	
	int cinfoF = 0;
	int cinfoM = 0;
	String errorstr = "";
	
	boolean chrx = false;
	
	
	
	
	public NF_NucFam(String id, ArrayList Data, ArrayList Loci, int ErrInterval, int Kprint, boolean x, String outfile)
	{
		ID = id;
		chrx = x;
		//Loci structure : chr \t rsnumber \t cM \t bp
		
		AllPositions = new int[Loci.size()];
		
		for(int i=0; i<Loci.size(); i++)
		{
			String[] posString = ((String)Loci.get(i)).split("\t");
			if(posString.length!=4) exit("Problem in list of loci : line length = "+posString.length+".");
			AllPositions[i] = Integer.parseInt(posString[3]);
		}
		
	
		//Data structure : ArrayList of Strings. Name+"-"Sex+"-"+SNP1+\t+SNP2+....+\t+SNPn et SNPx = allele1+" "+allele2
		//First line = father, Second line = mother and following lines = children
		
		//Parents
		
		String[] fatherData = ((String)Data.get(0)).split("-");
		
		if(fatherData.length != 3) exit("Problem in Data Format for Father in Nuclear Family "+ID);
		if(fatherData[2].split("\t").length != AllPositions.length) 
			exit("Problem in SNP Data Format for Father in Nuclear Family "+ID+" ("+fatherData[2].split("\t").length+" SNPs found): "+ AllPositions.length+" positions in map file");
		if(Integer.parseInt(fatherData[1])!=1) exit("Father is not coded as male (1) in Nuclear Family "+ID);
		
		Father = new NF_Parent(fatherData[0],fatherData[2], 1, this);
		
		
		String[] motherData = ((String)Data.get(1)).split("-");
		
		if(motherData.length != 3) exit("Problem in Data Format for Mother in Nuclear Family "+ID);
		if(motherData[2].split("\t").length != AllPositions.length) 
			exit("Problem in SNP Data Format for Mother in Nuclear Family "+ID+" ("+fatherData[2].split("\t").length+" SNPs found): "+ AllPositions.length+" positions in map file");
		if(Integer.parseInt(motherData[1])!=2) exit("Mother is not coded as female (2) in Nuclear Family "+ID);
		
		Mother = new NF_Parent(motherData[0],motherData[2], 2, this);
		
		locisize = Loci.size();
		//System.err.println("\n***** Data information for Family "+ ID +":");
		//System.err.println("* "+Loci.size()+" markers are found in data.");
		//System.err.println("* Loading parents ... done. ");
		
		
		SetHapInformative();
		//PrintHapInformative(outfile+"infomarkM.txt",'M');
		//PrintHapInformative(outfile+"infomarkF.txt",'F');
		
		
		//Children
		
		boolean suite = true;
		if(Data.size()<4){
			System.err.println("* Loading children ... less than 2 children in Family "+ID+" : impossible to compute crossover events!\n*****\n");
			suite = false;
		}
		
		
		
		if(suite)
		{
			Children = new TreeSet();
			
			int nbchild = 0;
			for(int i=2; i<Data.size(); i++)
			{
				nbchild++;
				String[] childData = ((String)Data.get(i)).split("-");
				if(childData.length != 3) exit("Problem in Data Format for Child "+nbchild+" in Nuclear Family "+ID);
				if(childData[2].split("\t").length != AllPositions.length) exit("Problem in SNP Data Format for Child "+nbchild+" in Nuclear Family "+ID);
				
				NF_Child act = new NF_Child(childData[0],childData[2],Integer.parseInt(childData[1]),(Data.size() - 2), this, chrx);
				Children.add(act);
			}
			
			//System.err.println("* Loading children ... "+Children.size()+" children found\n*****\n");
			
            //Print informative marker allele for each children
            
            PrintHapInformativeChildren(outfile+"infomarkM.txt",'M');
            PrintHapInformativeChildren(outfile+"infomarkF.txt",'F');
            
			//Recode each child with respect to the siblings + crossovers without checking for errors.
			
			/*
			System.err.println("Crossovers without looking for genotyping errors :\n");
			for (Iterator itr = Children.iterator(); itr.hasNext();)
			{
				NF_Child act = (NF_Child)itr.next();
				act.Recoding();
				act.Crossovers();
				
				
				//Print:
				
				act.printCrossovers('F');
				act.printCrossovers('M');
			}*/
			
			for (Iterator itr = Children.iterator(); itr.hasNext();)
			{
				NF_Child act = (NF_Child)itr.next();
				act.Recoding();
			}
			
			
			FindAllErrors(ErrInterval,false,outfile+"allErrors.txt");
			PrintErr(outfile+"allErrors.txt");
			
			//System.err.println("\nCROSSOVERS (after removing genotyping errors) :\n");
			
			Printing = new String[Kprint+2]; //+2 : legend + Kprint + Kprint=0
			Printing[0] = familyInfo()+"\nCROSSOVERS (after removing genotyping errors) :\n\n"; //always!
			for(int i=1; i<Kprint+2; i++)
			{
				Printing[i] = "";
			}
			
			for (Iterator itr = Children.iterator(); itr.hasNext();)
			{
				NF_Child act = (NF_Child)itr.next();
				
				act.Recoding();
				
				//Sharing
				if(Recomb.sbool)
				{
					String filesharing = outfile + "sharing.txt";
					
					ArrayList PaternalSharing = act.printRecodeErrorArray('F');
					ArrayList MaternalSharing = act.printRecodeErrorArray('M');
					PrintSharing(PaternalSharing, filesharing);
					PrintSharing(MaternalSharing, filesharing);
					
				}
				
				//Crossovers
				if(Recomb.cbool)
				{
				act.CrossoversErrors();
				
				//Print:
				
				//to have no Kprint
				
				for(int i=1; i<Kprint+2; i++)
				{
					//Printing[i] += act.printCrossovers('F', i) + act.printCrossovers('M', i);
					//to have no Kprint as 1
					if(!chrx) Printing[i] += act.printCrossovers('F', i-1) + act.printCrossovers('M', i-1);
					else Printing[i] += act.printCrossovers('M', i-1);
					
				}
				}
				
			}
			
		}
	
	}
	
	
	
	//Parent Management
	
	public String familyInfo()
	{
		String a = "";
		a += "\n***** Data information for Family "+ ID +":\n";
		a += "* "+locisize+" markers are found in data.\n";
		a += "* Loading parents ... done.\n";
		a += "* Father "+Father.Name+" presents "+cinfoF+" informative markers\n"; //("+Father.HapInformative+")");
		a += "* Mother "+Mother.Name+" presents "+cinfoM+" informative markers\n"; //("+Mother.HapInformative+")");
		a += "* Loading children ... "+Children.size()+" children found\n*****\n\n";
		a += errorstr;
		return a;
	}
	
	//add for cancer
	public void PrintErr(String file)
	{
		if(AllErrors.size()!=0)
		{
		try 
		{
			System.err.print("Writing in "+file+ " for family "+ID);
			BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
		
			bw.write("MarkerNbE\tPosition\tFam\n");
			for (Iterator itr = AllErrors.iterator(); itr.hasNext();)
			{
				Integer A = ((Integer)itr.next());
				bw.write(A+"\t"+AllPositions[A.intValue()]+"\t"+ID+"\n");
			}
			System.err.println(" done.");
			bw.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File writing problem");
		}
		}
	}
	
	//add for cancer
	//formatting changed for LVOTONULL part of cancer project
	
	public void PrintHapInformative(String file, char parent)
	{
		if(Father == null || Mother == null) exit("Impossible to use the PrintHapInformative method : parents are not defined (null)");
		
		if(parent == 'F')
		{
			boolean exit=false;
			if(Father.HapInformativeIsSet()==false)
			{
				System.err.println("Father's informative haplotype is not set");
				exit = true;
			}
			if(!exit)
			{
				if(chrx)System.err.println("No father's informative markers on X chromosome. No file is created.");
				else 
				{
					try 
					{
						System.err.print("Writing in "+file+ " for family "+ID);
						BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
						
						//bw.write("MarkerNbF\tPosition\tFam\n");
						for(int i=0; i<Father.HapInformative.length();i++)
						{
							if(Father.HapInformative.charAt(i)!='r')
							{
								bw.write(i+"\t"+AllPositions[i]+"\t"+ID+"\n");
                                
							}
						}
						bw.write("\nEnd\n");
						System.err.println(" done.");
						bw.close();
					}
					catch(Exception e)
					{
						e.printStackTrace();
						exit("File writing problem");
					}
					
				}

			}
		}
		else if(parent == 'M')
		{
			boolean exit=false;
			if(Mother.HapInformativeIsSet()==false)
			{
				System.err.println("Mother's informative haplotype is not set");
				exit = true;
			}
			if(!exit)
			{
				
				try 
				{
					System.err.print("Writing in "+file+ " for family "+ID);
					BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
					
					//bw.write("MarkerNbM\tPosition\tFam\n");
					for(int i=0; i<Mother.HapInformative.length();i++)
					{
						if(Mother.HapInformative.charAt(i)!='r')
						{
							bw.write(i+"\t"+AllPositions[i]+"\t"+ID+"\n");
						}
					}
					
					System.err.println(" done.");
					bw.close();
				}
				catch(Exception e)
				{
					e.printStackTrace();
					exit("File writing problem");
				}
								
			}
		}
		
	}
    
    //added for Imputation project
    
    public void PrintHapInformativeChildren(String file, char parent)
    {
        if(Father == null || Mother == null) exit("Impossible to use the PrintHapInformative method : parents are not defined (null)");
        
        
        
        
        
        if(parent == 'F')
        {
            boolean exit=false;
            if(Father.HapInformativeIsSet()==false)
            {
                System.err.println("Father's informative haplotype is not set");
                exit = true;
            }
            if(!exit)
            {
                if(chrx)System.err.println("No father's informative markers on X chromosome. No file is created.");
                else
                {
                    try
                    {
                        System.err.print("Writing in "+file+ " for family "+ID);
                        BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
                        
                        //bw.write("MarkerNbF\tPosition\tFam\n");
                        
                        bw.write("#ChildrensNames: ");
                        for (Iterator itr = Children.iterator(); itr.hasNext();)
                        {
                            NF_Child act = (NF_Child)itr.next();
                            bw.write(act.Name+" ");
                        }
                        bw.write("\n");
                        
                        for(int i=0; i<Father.HapInformative.length();i++)
                        {
                            if(Father.HapInformative.charAt(i)!='r')
                            {
                                bw.write(i+"\t"+AllPositions[i]+"\t"+ID+"\t");
                                
                                for (Iterator itr = Children.iterator(); itr.hasNext();)
                                {
                                    NF_Child act = (NF_Child)itr.next();
                                    bw.write(act.HapInfoFather.charAt(i)+" ");
                                }
                                bw.write("\n");
                            }
                        }
                        //bw.write("\nEnd\n");
                        System.err.println(" done.");
                        bw.close();
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                        exit("File writing problem");
                    }
                    
                }
                
            }
        }
        else if(parent == 'M')
        {
            boolean exit=false;
            if(Mother.HapInformativeIsSet()==false)
            {
                System.err.println("Mother's informative haplotype is not set");
                exit = true;
            }
            if(!exit)
            {
                
                try 
                {
                    System.err.print("Writing in "+file+ " for family "+ID);
                    BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
                    
                    //bw.write("MarkerNbM\tPosition\tFam\n");
                    bw.write("#ChildrensNames: ");
                    for (Iterator itr = Children.iterator(); itr.hasNext();)
                    {
                        NF_Child act = (NF_Child)itr.next();
                        bw.write(act.Name+" ");
                    }
                    bw.write("\n");
                    
                    for(int i=0; i<Mother.HapInformative.length();i++)
                    {
                        if(Mother.HapInformative.charAt(i)!='r')
                        {
                            bw.write(i+"\t"+AllPositions[i]+"\t"+ID+"\t");
                            
                            for (Iterator itr = Children.iterator(); itr.hasNext();)
                            {
                                NF_Child act = (NF_Child)itr.next();
                                bw.write(act.HapInfoMother.charAt(i)+" ");
                            }
                            bw.write("\n");
                        }
                    }
                    
                    System.err.println(" done.");
                    bw.close();
                }
                catch(Exception e)
                {
                    e.printStackTrace();
                    exit("File writing problem");
                }
                
            }
        }
        
    }
	
	//add for sharing 03-12-2013
	public void PrintSharing(ArrayList Sharing, String file)
	{
		try 
		{
			//System.err.print("Writing in "+file+ " for family "+ID);
			BufferedWriter bw = new BufferedWriter( new FileWriter( file, true ) );
			
			for (int i=0;i<Sharing.size();i++)
			{
				bw.write(Sharing.get(i)+"\n");
			}
			//System.err.println(" done.");
			bw.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			exit("File writing problem in PrintSharing()");
		}
		
	}
	
	public void SetHapInformative()
	{
		if(Father == null || Mother == null) exit("Impossible to use the SetHapInformative method : parents are not defined (null)");
		
		boolean exit = false;
		if(Father.HapInformativeIsSet())
		{
			System.err.println("Father's informative haplotype is already set");
			exit = true;
		}
		if(Mother.HapInformativeIsSet())
		{
			System.err.println("Mother's informative haplotype is already set");
			exit = true;
		}
		
		
		if(!exit)
		{
			int countinfoF = 0;
			int countinfoM = 0;
			int missingF = 0;
			int missingM = 0;
			
			for(int i=0; i<Father.Data.length;i++)
			{
				String[] nucF = Father.Data[i].split(" ");
				String[] nucM = Mother.Data[i].split(" ");
				
				if(nucF.length!=2 || nucM.length!=2) exit("Problem in SNP Data Format (required : all1+\" \"+all2)");
				
				if((nucM[0].charAt(0)=='0' && nucM[1].charAt(0)=='0') || (nucF[0].charAt(0)=='0' && nucF[1].charAt(0)=='0'))
				{
					if(!chrx) Father.HapInformative += 'r';
					Mother.HapInformative += 'r';
				}
				
				else{
					
				if(!chrx)
				{
				//Father informative markers : Father is heterozygous and Mother is homozygous
				
				if(nucM[0].charAt(0)==nucM[1].charAt(0) && nucF[0].charAt(0)!=nucF[1].charAt(0)) //M homoz F heteroz
				{
					//no missing data in parents
					if(nucM[0].charAt(0)!='0' && nucM[1].charAt(0)!='0' && nucF[0].charAt(0)!='0' && nucF[1].charAt(0)!='0')
					{
						countinfoF++;
						if(nucM[0].charAt(0)==nucF[0].charAt(0))
						{
							Father.HapInformative += nucF[1].charAt(0);
						}
						else if(nucM[0].charAt(0)==nucF[1].charAt(0)) //misssing
						{
							Father.HapInformative += nucF[0].charAt(0);
						}
					}
					//missing data in parents = not informative
					else
					{
						//missingM++;
						Father.HapInformative += 'r';
					}
				}
				else //not informative
				{
					Father.HapInformative += 'r';
				}
				
				}
				
				
				
				//Mother informative markers : Mother is heterozygous and Father is homozygous
				
				
				if(nucM[0].charAt(0)!=nucM[1].charAt(0) && nucF[0].charAt(0)==nucF[1].charAt(0))
				{
					if(nucM[0].charAt(0)!='0' && nucM[1].charAt(0)!='0' && nucF[0].charAt(0)!='0' && nucF[1].charAt(0)!='0')
					{
						countinfoM++;
						if(nucF[0].charAt(0)==nucM[0].charAt(0))
						{
							Mother.HapInformative += nucM[1].charAt(0);
						}
						else if(nucF[0].charAt(0)==nucM[1].charAt(0)) //misssing
						{
							Mother.HapInformative += nucM[0].charAt(0);
						}
					}
					else
					{
						//missingF++;
						Mother.HapInformative += 'r';
					}
				}
				else
				{
					Mother.HapInformative += 'r';
				}
				
				}
				
			}
			
			cinfoF = countinfoF;
			cinfoM = countinfoM;
			
			if(chrx) Father.HapInformative = "-";
			
			//System.err.println("* Father "+Father.Name+" presents "+countinfoF+" informative markers"); //("+Father.HapInformative+")");
			//System.err.println("* Mother "+Mother.Name+" presents "+countinfoM+" informative markers"); //("+Mother.HapInformative+")");
		}
		
	}
	
	//Recup of all errors from children in TreeSet AllErrors.
	
	public void FindAllErrors(int ErrInterval, boolean printall, String file)
	{
		
		//System.err.println("ERRORS (if ErrInterval is < "+ ErrInterval +") :\n");
		errorstr = "ERRORS (if ErrInterval is < "+ ErrInterval +") written in "+file+"\n";
		if(printall) errorstr +="\n";

		
		AllErrors = new TreeSet();
		
		for (Iterator itr = Children.iterator(); itr.hasNext();)
		{
			NF_Child act = (NF_Child)itr.next();
			if (act.recode)
			{
			
				String res = act.FindErrors(ErrInterval); //act.Errors is defined.
				if(printall)errorstr += res;
				for (Iterator itr2 = act.Errors.iterator(); itr2.hasNext();)
				{
					Integer A = ((Integer)itr2.next());
					AllErrors.add(A);
				}
			
			}
			else
			{
				exit("Cannot use FindAllErrors() : at least one child is not recoded");
			}
		}
		
		errorstr += "\n\tTotal:\t"+AllErrors.size()+" errors found in children for Family "+ ID +".\n";
	
	}
	
	
	
	
	
	//Utilitaires
	
	public static void exit(String s)
	{
	
		System.err.println(s);
		System.exit(0);
	
	}
	


}
