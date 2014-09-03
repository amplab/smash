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

import java.util.function.Function;

/**
 * A type used for submarining checked exceptions through methods that don't declare them, for
 * example, {@link Function#apply}. The top level {@link Main#main} method contains a
 * {@code try-catch} block to catch and unwrap these exceptions.
 */
public class ExceptionWrapper extends RuntimeException {

  public static ExceptionWrapper wrap(Exception cause) {
    return new ExceptionWrapper(cause);
  }

  private ExceptionWrapper(Exception cause) {
    super(cause);
  }

  @Override
  public Exception getCause() {
    return (Exception) super.getCause();
  }
}
