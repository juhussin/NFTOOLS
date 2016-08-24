//
//  NF_Parent.java
//  LVTOOLS
//
//  Created by Julie Hussin on 09-07-31.
//  Copyright 2009 __MyCompanyName__. All rights reserved.

//
//	Object NF_Parent
//

//package NUCFAMTOOLS;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class NF_Parent {

	//Identity
	String Name;
	int Sex;			//1=male, 2=female
	NF_NucFam FamilyP;	//this is a parent in FamilyP
	NF_NucFam FamilyC;	//this is a child in FamilyC
	
	//Data
	String[] Data;
	String HapInformative;	//Haplotype for informative markers with respect to NF_Parent FamilyP.otherParent
							//Nucleotide when informative, 'r' when not.
							
							
	public NF_Parent(String name, String data, int sex, NF_NucFam famP)
	{
	
		//Data structure : SNP1+\t+SNP2+....+\t+SNPn
		Name = name;
		Data = data.split("\t");
		Sex = sex;	
		FamilyP = famP;
		FamilyC = null;
		HapInformative = "";	//has to be set via NucFam class
	
	}
	
	public boolean HapInformativeIsSet(){
	
		return !HapInformative.equals("");
	
	}
	
	public void SetFamilyC(NF_NucFam family){
	
		FamilyC = family;
	
	}
	
	public void SetHapInformative(){
	
		FamilyP.SetHapInformative();
	
	}
	

}