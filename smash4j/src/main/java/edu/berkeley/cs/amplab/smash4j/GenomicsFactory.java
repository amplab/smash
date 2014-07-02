package edu.berkeley.cs.amplab.smash4j;

import com.google.api.client.extensions.java6.auth.oauth2.AuthorizationCodeInstalledApp;
import com.google.api.client.extensions.jetty.auth.oauth2.LocalServerReceiver;
import com.google.api.client.googleapis.auth.oauth2.GoogleAuthorizationCodeFlow;
import com.google.api.client.googleapis.auth.oauth2.GoogleClientSecrets;
import com.google.api.client.googleapis.auth.oauth2.GoogleCredential;
import com.google.api.client.googleapis.javanet.GoogleNetHttpTransport;
import com.google.api.client.googleapis.services.CommonGoogleClientRequestInitializer;
import com.google.api.client.googleapis.services.GoogleClientRequestInitializer;
import com.google.api.client.http.HttpRequest;
import com.google.api.client.http.HttpRequestInitializer;
import com.google.api.client.http.HttpTransport;
import com.google.api.client.json.JsonFactory;
import com.google.api.client.json.jackson2.JacksonFactory;
import com.google.api.client.util.store.DataStoreFactory;
import com.google.api.client.util.store.FileDataStoreFactory;
import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.GenomicsScopes;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.security.GeneralSecurityException;
import java.util.Collection;
import java.util.prefs.Preferences;

/**
 * A factory for {@link Genomics} stubs.
 */
public final class GenomicsFactory {

  /**
   * A builder for {@link GenomicsFactory} objects.
   */
  public final static class Builder {

    private String applicationName;
    private DataStoreFactory dataStoreFactory;
    private HttpTransport httpTransport;
    private JsonFactory jsonFactory;
    private Collection<String> scopes;
    private String userName;
    private int readTimeout;
    private int connectTimeout;

    private Builder() throws GeneralSecurityException, IOException {
      setApplicationName(Main.PROGRAM_NAME);
      setDataStoreFactory(
          new FileDataStoreFactory(new File(System.getProperty("user.home"), ".store/smash4j")));
      setHttpTransport(GoogleNetHttpTransport.newTrustedTransport());
      setJsonFactory(JacksonFactory.getDefaultInstance());
      setScopes(GenomicsScopes.all());
      setUserName(System.getProperty("user.name"));
      setReadTimeout(20000);
      setConnectTimeout(20000);
    }

    /**
     * Set the application name. Default is {@code AMPLab-SMaSH4J/0.1}.
     */
    public Builder setApplicationName(String applicationName) {
      this.applicationName = applicationName;
      return this;
    }

    /**
     * Set the {@link DataStoreFactory}. Default is {@code ~/.store/smash4j}.
     */
    public Builder setDataStoreFactory(DataStoreFactory dataStoreFactory) {
      this.dataStoreFactory = dataStoreFactory;
      return this;
    }

    /**
     * Set the {@link HttpTransport}. Default is {@link GoogleNetHttpTransport#newTrustedTransport}.
     */
    public Builder setHttpTransport(HttpTransport httpTransport) {
      this.httpTransport = httpTransport;
      return this;
    }

    /**
     * Set the {@link JsonFactory}. Default is {@link JacksonFactory#getDefaultInstance}.
     */
    public Builder setJsonFactory(JsonFactory jsonFactory) {
      this.jsonFactory = jsonFactory;
      return this;
    }

    /**
     * Set the scopes. Default is {@link GenomicsScopes#all}.
     */
    public Builder setScopes(Collection<String> scopes) {
      this.scopes = scopes;
      return this;
    }

    /**
     * Set the user name. Default is the system property {@code user.name}.
     */
    public Builder setUserName(String userName) {
      this.userName = userName;
      return this;
    }

    /**
     * Set the read timeout.
     */
    public Builder setReadTimeout(int readTimeout) {
      this.readTimeout = readTimeout;
      return this;
    }

    /**
     * Set the connect timeout.
     */
    public Builder setConnectTimeout(int connectTimeout) {
      this.connectTimeout = connectTimeout;
      return this;
    }

    /**
     * Build the {@link GenomicsFactory}.
     */
    public GenomicsFactory build() {
      return new GenomicsFactory(
          applicationName,
          dataStoreFactory,
          httpTransport,
          jsonFactory,
          scopes,
          userName,
          readTimeout,
          connectTimeout);
    }
  }

  /**
   * Static factory method for {@link Builder}.
   */
  public static Builder builder() throws GeneralSecurityException, IOException {
    return new Builder();
  }

  private static IllegalStateException badFile(File file, String key, String reason) {
    return illegalStateException(
        "File \"%s\" specified in %s:%s must %s.",
        file.getPath(), Main.PREFERENCES_PATH, key, reason);
  }

  private static File checkFile(String path, String key) {
    File file = new File(path);
    if (!file.exists()) {
      throw badFile(file, key, "exist");
    }
    if (!file.isFile()) {
      throw badFile(file, key, "be a file");
    }
    if (!file.canRead()) {
      throw badFile(file, key, "be readable");
    }
    return file;
  }

  /**
   * Return the default {@link GenomicsFactory}.
   */
  public static GenomicsFactory getDefault() throws GeneralSecurityException, IOException {
    return builder().build();
  }

  private static String getPreference(String key, String format, Object... args) {
    String value =  Main.PREFERENCES.get(key, null);
    if (null == value) {
      throw illegalStateException(String.format(format, args));
    }
    return value;
  }

  private static IllegalStateException illegalStateException(String format, Object... args) {
    return new IllegalStateException(String.format(format, args));
  }

  private final String applicationName;
  private final DataStoreFactory dataStoreFactory;
  private final HttpTransport httpTransport;
  private final JsonFactory jsonFactory;
  private final Collection<String> scopes;
  private final String userName;
  private final int readTimeout;
  private final int connectTimeout;

  private GenomicsFactory(
      String applicationName,
      DataStoreFactory dataStoreFactory,
      HttpTransport httpTransport,
      JsonFactory jsonFactory,
      Collection<String> scopes,
      String userName,
      int readTimeout,
      int connectTimeout) {
    this.applicationName = applicationName;
    this.dataStoreFactory = dataStoreFactory;
    this.httpTransport = httpTransport;
    this.jsonFactory = jsonFactory;
    this.scopes = scopes;
    this.userName = userName;
    this.readTimeout = readTimeout;
    this.connectTimeout = connectTimeout;
  }

  /**
   * Create a new {@link Genomics} stub from the given API key.
   */
  public Genomics fromApiKey(String apiKey) {
    return create(null, new CommonGoogleClientRequestInitializer(apiKey));
  }

  /**
   * Create a new {@link Genomics} stub from the given client secrets JSON file.
   */
  public Genomics fromClientSecretsFile(File clientSecretsJson) throws IOException {
    try (Reader in = new FileReader(clientSecretsJson)) {
      GoogleAuthorizationCodeFlow flow = new GoogleAuthorizationCodeFlow
          .Builder(
              httpTransport,
              jsonFactory,
              GoogleClientSecrets.load(jsonFactory, in),
              scopes)
          .setDataStoreFactory(dataStoreFactory)
          .build();
      return create(
          new AuthorizationCodeInstalledApp(flow, new LocalServerReceiver()).authorize(userName),
          null);
    }
  }

  /**
   * Create a new {@link Genomics} stub from the values in {@link Preferences}.
   */
  public Genomics fromPreferences() throws GeneralSecurityException, IOException {
    String authorizationMethod = getPreference(
        "authorizationMethod",
        "%s:authorizationMethod is unset. Please set it to one of \"API_KEY\", "
            + "\"SERVICE_ACCOUNT\", or \"CLIENT_SECRETS\".",
        Main.PREFERENCES_PATH);
    switch (authorizationMethod) {
      case "API_KEY":
        return fromApiKey(
            getPreference(
                "apiKey",
                "When %s:authorizationMethod = \"API_KEY\", %s:apiKey must be set.",
                Main.PREFERENCES_PATH,
                Main.PREFERENCES_PATH));
      case "SERVICE_ACCOUNT":
        return fromServiceAccount(
            getPreference(
                "serviceAccountId",
                "When %s:authorizationMethod = \"SERVICE_ACCOUNT\", %s:serviceAccountId must be "
                    + "set.",
                Main.PREFERENCES_PATH,
                Main.PREFERENCES_PATH),
            checkFile(
                getPreference(
                    "serviceAccountP12File",
                    "When %s:authorizationMethod = \"SERVICE_ACCOUNT\", %s:serviceAccountP12File "
                        + "must be set.",
                    Main.PREFERENCES_PATH,
                    Main.PREFERENCES_PATH),
                "serviceAccountP12File"));
      case "CLIENT_SECRETS":
        return fromClientSecretsFile(
            checkFile(
                getPreference(
                    "clientSecretsFile",
                    "When %s:authorizationMethod = \"CLIENT_SECRETS\", %s:clientSecretsFile must "
                        + "be set.",
                    Main.PREFERENCES_PATH,
                    Main.PREFERENCES_PATH),
                "clientSecretsFile"));
      default:
        throw illegalStateException(
            "%s:authorizationMethod set to an invalid value \"%s\". Please set it to one of "
                + "{ \"API_KEY\", \"SERVICE_ACCOUNT\", \"CLIENT_SECRETS\" }.",
            Main.PREFERENCES_PATH,
            authorizationMethod);
    }
  }

  /**
   * Create a new {@link Genomics} stub from the given service account ID and P12 file.
   */
  public Genomics fromServiceAccount(String serviceAccountId, File p12File)
      throws GeneralSecurityException, IOException {
    return create(
        new GoogleCredential.Builder()
            .setTransport(httpTransport)
            .setJsonFactory(jsonFactory)
            .setServiceAccountId(serviceAccountId)
            .setServiceAccountScopes(scopes)
            .setServiceAccountPrivateKeyFromP12File(p12File)
            .build(),
        null);
  }

  private Genomics create(
      final HttpRequestInitializer delegate,
      GoogleClientRequestInitializer googleClientRequestInitializer) {
    return new Genomics
        .Builder(
            httpTransport,
            jsonFactory,
            new HttpRequestInitializer() {
              @Override public void initialize(HttpRequest httpRequest) throws IOException {
                if (null != delegate) {
                  delegate.initialize(httpRequest);
                }
                httpRequest.setReadTimeout(readTimeout);
                httpRequest.setConnectTimeout(connectTimeout);
              }
            })
        .setApplicationName(applicationName)
        .setGoogleClientRequestInitializer(googleClientRequestInitializer)
        .build();
  }
}
