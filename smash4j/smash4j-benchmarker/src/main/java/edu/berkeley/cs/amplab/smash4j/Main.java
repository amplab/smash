package edu.berkeley.cs.amplab.smash4j;

import com.google.api.services.genomics.Genomics;

import java.io.IOException;
import java.security.GeneralSecurityException;

public class Main {

  public static void main(String[] args) throws GeneralSecurityException, IOException {
    Genomics genomics = GenomicsFactory.getDefault().fromPreferences();
  }
}
