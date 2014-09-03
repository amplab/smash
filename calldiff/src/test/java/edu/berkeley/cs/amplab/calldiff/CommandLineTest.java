/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package edu.berkeley.cs.amplab.calldiff;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import java.util.Optional;

import edu.berkeley.cs.amplab.calldiff.CommandLine;

/**
 * Unit test for {@link CommandLine}
 */
public class CommandLineTest {

  @Test
  public void testCommandLine() {
    assertEquals(
        Optional.of(CommandLine.builder()
            .setLhsSampleId("lhs_sample_id")
            .setLhsVcf("lhs_vcf")
            .setPresorted(true)
            .setReferenceFai("reference_fai")
            .setReferenceFasta("reference_fasta")
            .setRhsSampleId("rhs_sample_id")
            .setRhsVcf("rhs_vcf")
            .build()),
        CommandLine.parse(
            "--lhs_sample_id=lhs_sample_id",
            "--lhs_vcf=lhs_vcf",
            "--presorted",
            "--reference_fai=reference_fai",
            "--reference_fasta=reference_fasta",
            "--rhs_sample_id=rhs_sample_id",
            "--rhs_vcf=rhs_vcf"));
  }
}
