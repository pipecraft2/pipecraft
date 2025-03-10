require('dotenv').config();
const { notarize } = require('@electron/notarize');

exports.default = async function notarizing(context) {
  const { electronPlatformName, appOutDir } = context;
  if (electronPlatformName !== 'darwin') {
    return;
  }

  // Check if environment variables are set
  if (!process.env.APPLE_ID || !process.env.APPLE_APP_SPECIFIC_PASSWORD) {
    throw new Error('APPLE_ID and APPLE_APP_SPECIFIC_PASSWORD must be set in .env file');
  }

  const appName = context.packager.appInfo.productFilename;

  console.log('Notarizing application...');
  console.log(`App path: ${appOutDir}/${appName}.app`);
  
  try {
    await notarize({
      appBundleId: 'com.pipecraft.app',
      appPath: `${appOutDir}/${appName}.app`,
      appleId: process.env.APPLE_ID,
      appleIdPassword: process.env.APPLE_APP_SPECIFIC_PASSWORD,
      teamId: '3CYD3SH29Q'
    });
  } catch (error) {
    console.error('Notarization failed:', error);
    throw error;
  }
  
  console.log('Notarization complete!');
}; 